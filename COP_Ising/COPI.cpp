/* COMPILE WITH:
   g++ -fopenmp -O3 -march=native -std=c++20 COPI.cpp -o COPI

   RUN WITH:
   ./COPI [int]
*/

/************************************************************************************************
  This program simulates the N-dimensional Ising model out of equilibrium with the enhanced bulk 
  diffusion algorithm (Marko and Barkema, 1995) and generalized helical boundary conditions.
 ************************************************************************************************/

//NOTE: To use periodic boundary conditions, make S=0. To use ordinary screw-periodic (or helical) boundary conditions, make S=1. 
// If you wish to use generalized screw-periodic boundary conditions, set S to any integer between 2 and L-1 (L is the lattice size).  
// If (S mod L = 0), then this is just periodic boundary conditions, else it's screw-periodic boundary conditions with S = (S mod L).

#include <iostream>   // Input and output
#include <vector>     // Dynamic arrays
#include <cmath>      // Sin, cos, exp, pow, sqrt, atan2, ...
#include <fstream>    // Reading and writing to files
#include <string>     // Unique file names 
#include <iomanip>    // Set precision of output data
#include <algorithm>  // Std::fill (a little faster than for loop)
#include <cstdint>    // Fixed width integer types (in prng and Wolff)
#include <chrono>     // Getting current time
#include <ctime>      // Making time from chrono more human readable
#include <sstream>    // Treats strings as files 
#include <omp.h>      // Parallelize
#include <filesystem> // To check if file exists
#include <array>      // to get fixed array sizes at compile time (makes program run roughly 40% faster, but compilation is a little slower)

using RealType = double;
using IntType = long long;

/* xoshiro256+ implementation by David Blackman and Sebastiano Vigna (vigna@acm.org)
 * Public domain.
 *
 * The period of this prng is 2^256 - 1. 
 */
struct Xoshiro256Plus { 
    uint64_t state[4];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    // splitmix64 for seeding as recommended by the authors of xoshiro256+
    void seed(uint64_t seed_val) {
        uint64_t z = seed_val;
        for (int i = 0; i < 4; i++) {
            z += 0x9e3779b97f4a7c15;
            uint64_t s = z;
            s = (s ^ (s >> 30)) * 0xbf58476d1ce4e5b9;
            s = (s ^ (s >> 27)) * 0x94d049bb133111eb;
            state[i] = s ^ (s >> 31);
        }
    }

    uint64_t next() {
        const uint64_t result = state[0] + state[3];
        const uint64_t t = state[1] << 17;
        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];
        state[2] ^= t;
        state[3] = rotl(state[3], 45);
        return result;
    }

    double next_double() {
        return (double)(next() >> 11) * 0x1.0p-53;
    }

    // Advances 2^128 steps.
    void jump() {
        static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < 4; i++) {
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & 1ULL << b) {
                    s0 ^= state[0];
                    s1 ^= state[1];
                    s2 ^= state[2];
                    s3 ^= state[3];
                }
                next();
            }
        }
        state[0] = s0;
        state[1] = s1;
        state[2] = s2;
        state[3] = s3;
    }

    // Advances 2^192 steps.
    void long_jump() {
        static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe69d };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < 4; i++) {
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & 1ULL << b) {
                    s0 ^= state[0];
                    s1 ^= state[1];
                    s2 ^= state[2];
                    s3 ^= state[3];
                }
                next();
            }
        }
        state[0] = s0;
        state[1] = s1;
        state[2] = s2;
        state[3] = s3;
    }
};

inline int fast_mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

// represents an instance of the system being simulated
template <int dim>   // to speed up simulation for common instances of dim (e.g., dim=2,3 predominantly get used as opposed to dim=912 for instance)
struct alignas(64) System { // alignas(64) guarantees no false sharing occurs later on
    static constexpr int z = 2 * dim + 1; // lattice coordination number
    int S;      // parameter for generalized helical boundary conditions
    int L;      // number of sites along a single dimension (assumed equal)
    int Total;  // Total number of sites 
    RealType t; // total time elapsed for this system
    RealType beta; // inverse temperature

    std::vector<int> Neighbors;      // data structure to access neighbors of every site. 
                                     // First entry is in positive direction of dimension N, 
                                     // second entry is in positive direction of dimension N-1, ... ,
                                     // the (N+1)th entry is in the negative direction of dimension N, 
                                     // the (N+2)th entry is in the negative direction of dimension N-1, ...
    std::vector<int8_t> Status;      // vector that stores spin configuration of every site (+/- 1)
    std::vector<int> neighb_Q;       // stores information where a certain site address is in the Q vector
    std::vector<int> Q;              // Q[i*Total + j] is the address of the j-th spin with coordination number i 

    // fixed size at compile time (these depent only on dim)
    std::array<int, 2*dim> neighb_flip;        // in Iterate(), this stores neighbors of a site that are antialigned with it
    std::array<int, dim> positions;            // auxiliary to set up boundary conditions (stores coordinates)
    std::array<int, dim> temp_positions;       // another auxiliary to set up boundary conditions 
    std::array<int, z> NQ;                     // NQ[i] is the total number of spins with coordination number i
    std::array<int, 4*dim - 2> update_neighb;  // stores all neighbors that will have coordinatioon number changed after a flip
    std::array<RealType, z> factor;            // precalculated term (1 - q/z)*Nq*exp(-4Jq) in Marko and Barkema's paper 
    Xoshiro256Plus rng;                        // Each system has its own PRNG

    void init_sys(int L_in, int S_in, bool random_start, RealType beta_in) { // (assumes prng is seeded)
        L     = L_in;
        S     = S_in;
        t     = 0.0;
        beta  = beta_in;
        Total = static_cast<int>(std::pow(L, dim));
        Neighbors.resize(Total * 2 * dim); // only nearest-neighbors
        Status.resize(Total);              
        Q.resize(Total*z);     // maximum in theory is every site has same coordination number 
        neighb_Q.resize(Total);

        // Initialize spins
        if (random_start) { // hot start
            for (int i = 0; i < Total; ++i) {
                Status[i] = 2*(static_cast<int>(2.0*rng.next_double())) - 1; // -1 or +1 with equal probability
            }
        } else { // cold start 
            for (int i = 0; i < Total; ++i) { 
                Status[i] = 1; 
            }
        }

        // Initialize other data structures
        for (int i = 0; i < z; ++i) {
            factor[i] = (1.0 - static_cast<RealType>(i) / static_cast<RealType>(2*dim))*std::exp(-4.0*i*beta);
            NQ[i]     = 0;
        }

        set_neighbors_recursive(0);
        for (int i = 0; i < Total; ++i){  // for every site... 
            int coord_num = 0; 
            for (int j = 0; j < 2*dim; ++j) { // ... for each of its neighbors... 
                if (Status[i] + Status[Neighbors[2*dim*i + j]]){ // ... if they're aligned, coordination number grows
                    ++coord_num;                                 // (if statement evaluates to True if argument is anything but zero)
                }
            }
            neighb_Q[i]    = coord_num*Total + NQ[coord_num];
            Q[neighb_Q[i]] = i; 
            ++(NQ[coord_num]); // acts sort of like a stack pointer
        }
    }

    // maps n dimensional vector: (x_1, x_2, ..., x_dim) , where x_i = 0, 1, ..., L-1 to positive integers. 
    // e.g., (0,0) -> 0 , (0,1) -> 1 , ... , (1,0) -> L , ... , (L-1,L-1) -> L^2 - 1
    int f() const {
        int res = 0; 
        for (int i = 0; i < dim; ++i) res += positions[i] * static_cast<int>(std::pow(L, dim - 1 - i)); 
        return res;
    }

    void set_neighbors_recursive(int depth) { 
        if (depth >= dim) { // found a viable combination
            int temp = f();
            for (int j = 0; j < dim; ++j) {
                for (int k = 0; k < dim; ++k) 
                    temp_positions[k] = positions[k]; // save state 
                
                // "forward" direction in that dimension
                int x = (positions[j] + 1) % L; 
                if (x != 0) { // not on boundary
                    positions[j] = x; 
                    Neighbors[temp * (2*dim) + j] = f(); 
                } else { // on boundary
                    positions[j] = 0; 
                    int tempj    = j;
                    for (int k = 0; k < (dim - 1); ++k) { // for other coordinates in cyclic order
                        tempj = (j + k + 1) % dim; 
                        x     = (positions[tempj] + S) % L; // S = 0 is PBC, S = 1 is normal SBC
                        if (x != 0) {
                            positions[tempj] = x; 
                            break;
                        } else {
                            positions[tempj] = 0; 
                        }
                    }
                    Neighbors[temp * (2*dim) + j] = f();
                }
                for (int k = 0; k < dim; ++k) 
                    positions[k] = temp_positions[k]; // put back state

                // "backward" direction in that dimension
                x = fast_mod(positions[j] - 1, L);
                if (x != (L - 1)) { // not on boundary
                    positions[j] = x; 
                    Neighbors[temp * (2*dim) + (j + dim)] = f(); 
                } else { // on boundary
                    positions[j] = L - 1; 
                    int tempj    = j;
                    for (int k = 0; k < (dim - 1); ++k) { // for other coordinates in cyclic order
                        tempj = (j + k + 1) % dim; 
                        x     = fast_mod(positions[tempj] - S, L); 
                        if (x != (L - 1)) {
                            positions[tempj] = x; 
                            break;
                        } else {
                            positions[tempj] = L - 1; 
                        }
                    }
                    Neighbors[temp * (2*dim) + (j + dim)] = f();
                }
                for (int k = 0; k < dim; ++k) 
                    positions[k] = temp_positions[k]; // put back state
            }
        } else { 
            for (int i = 0; i < L; ++i) { // generates all combinations of (x_1, x_2, ..., x_dim), where x_i = 0, 1, 2, ..., L-1
                positions[depth] = i;
                set_neighbors_recursive(depth + 1);
            }
        }
    }
  
    // Helper function to remove an abstract site from its current bin structure
    inline void remove_from_bin(int target_site) {
        int target_abs_pos = neighb_Q[target_site];
        int target_bin     = target_abs_pos / Total;

        --(NQ[target_bin]);
        int replacement_site = Q[target_bin * Total + NQ[target_bin]];
        Q[target_abs_pos] = replacement_site;
        neighb_Q[replacement_site] = target_abs_pos;
    }

    // Helper function to insert an abstract site into its correct bin structure
    inline void insert_into_bin(int target_site, int target_bin) {
        int new_abs_pos = target_bin * Total + NQ[target_bin];
        Q[new_abs_pos]  = target_site;
        neighb_Q[target_site] = new_abs_pos;
        ++(NQ[target_bin]);
    }

    // Performs one iteration of the bulk diffusion algorithm
    void Iterate() {
        RealType dt = 0.0; 
        for (int i = 0; i < z; ++i) {
            dt += factor[i]*NQ[i];
        }
        dt = 1.0 / dt;
        t += dt;

        RealType r = rng.next_double();
        RealType prob = dt * NQ[0] * factor[0];
        int bin = 0; 
        
        // We only check up to z-2 because factor[z-1] (which is factor[2*dim]) is strictly 0.0
        while (bin < z - 2 && r > prob) {
            ++bin;
            prob += dt * NQ[bin] * factor[bin]; 
        }

        // randomly select a site from the chosen list
        int site_idx = int(rng.next_double() * NQ[bin]); 
        int site     = Q[bin*Total + site_idx]; 

        // randomly select an antialigned neighbor
        int antialigned_neighbors = 0;
        for (int neighb = 0; neighb < 2*dim; ++neighb){
            if (Status[site] - Status[Neighbors[site*2*dim + neighb]]) { 
                neighb_flip[antialigned_neighbors++] = Neighbors[site*2*dim + neighb]; 
            }
        }
        int neighbor = neighb_flip[static_cast<int>(antialigned_neighbors*rng.next_double())]; 

        // Now update the 2 spins that will be flipped and all their neighbors
        // first, since we have the 2 sites which will be flipped, gather all the other neighbors
        int np = 0;
        for (int neighb = 0; neighb < 2*dim; ++neighb) {
            int n1 = Neighbors[site*2*dim + neighb];
            if (n1 != neighbor){
                update_neighb[np++] = n1;
            }

            int n2 = Neighbors[neighbor*2*dim + neighb];
            if (n2 != site) {
                update_neighb[np++] = n2;
            }
        }
    
        // now remove these sites from the data structures 
        remove_from_bin(site); 
        remove_from_bin(neighbor);
        for (int i = 0; i < np; ++i){
            remove_from_bin(update_neighb[i]);
        }
        
        // flip the spins
        Status[site] = -Status[site];
        Status[neighbor] = -Status[neighbor];

        // re-insert the 2 flipped spins into data structures 
        int coord_num1 = 0; 
        int coord_num2 = 0; 
        for (int neighb = 0; neighb < 2*dim; ++neighb) {
            if (Status[site] + Status[Neighbors[2*dim*site + neighb]]){
                ++coord_num1;
            }
            if (Status[neighbor] + Status[Neighbors[2*dim*neighbor + neighb]]){
                ++coord_num2;
            }
        }
        insert_into_bin(site, coord_num1);
        insert_into_bin(neighbor, coord_num2);
 
        // re-insert remaining affected neighbors too...
        for (int i = 0; i < np; ++i){
            int new_coord = 0;
            int affected_site = update_neighb[i];
            for (int neighb = 0; neighb < 2 * dim; ++neighb) {
                if (Status[affected_site] + Status[Neighbors[affected_site * 2 * dim + neighb]]) {
                    ++new_coord;
                }
            }
            insert_into_bin(affected_site, new_coord);
        }
    }

    // stores spin configuration at current time in simulation in specified output file
    void Save_config(const std::string& output_file_name) {
        std::ofstream out_append(output_file_name, std::ios::app);
        
        // Save simulation time, followed by the spins
        out_append << std::scientific << std::setprecision(8) << t << " ";
        for (int i = 0; i < Total; ++i) {
            // the + in front of status implicitely promotes int8_t to int. 
            // int8_t is understood as a signed char, so without this, the output would be wacky characters instaed of +1 and -1
            out_append << +Status[i] << (i == Total - 1 ? "" : " ");
        }
        out_append << "\n";
    }
};

/*
 * Runs an ensemble of independent COP Ising models with the bulk diffusion algorithm,
 * saving spin configurations at logarithmic time intervals.
 */
template <int dim>
void Run(int int_param, int seed, int L, RealType beta, int S, bool random_start, int M, RealType t_max, RealType log_base) { 
    int Nthreads = omp_get_max_threads();
    std::vector<System<dim>> systems(Nthreads); // template initialization
    
    // Master PRNG to seed each sample uniquely without overlap
    Xoshiro256Plus master_rng;
    master_rng.seed(seed + int_param);
    
    // Format beta to 3 decimal places for the filenames
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << beta;
    std::string beta_str = ss.str();
    
    // Pre-generate filenames
    std::vector<std::string> filenames(M);
    for (int i = 0; i < M; ++i) {
        filenames[i] = "configs_" + std::to_string(L) + "_" + beta_str + "_" + std::to_string(i) + ".dat";
    }

    int samples_simulated = 0;
    
    // Loop until all M samples are simulated
    while (samples_simulated < M) {
        // Only spin up threads for the remaining samples
        int current_batch_size = std::min(Nthreads, M - samples_simulated);
        
        // Initialize systems and save initial configurations
        for (int i = 0; i < current_batch_size; ++i) {
            int sample_idx = samples_simulated + i;
            
            // Copy master PRNG state to the system, then advance master
            for (int k = 0; k < 4; ++k) {
                systems[i].rng.state[k] = master_rng.state[k];
            }
            master_rng.jump(); 
            
            systems[i].init_sys(L, S, random_start, beta);
            
            // If file doesn't already exist, initialize it and save t=0.0 config
            if (!std::filesystem::exists(filenames[sample_idx])) {
                systems[i].Save_config(filenames[sample_idx]);
            }
        }
        
        // schedule(dynamic, 1) balances the load if different threads hit the save blocks at slightly different times
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < current_batch_size; ++i) {
            int sample_idx = samples_simulated + i;
            
            // First non-zero time to save.
            RealType next_save_t = 1.0; 
            
            while (systems[i].t < t_max) {
                systems[i].Iterate();
               
                // Check if we hit the next logarithmic interval
                if (systems[i].t >= next_save_t) {
                    #pragma omp critical  // once need to save, each thread will print current progress one by one (avoids mangled output)
                    {
                        std::cout << systems[i].t << std::endl;
                    }

                    systems[i].Save_config(filenames[sample_idx]);
                    next_save_t *= log_base; 
                }
            }
        }
        samples_simulated += current_batch_size;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) { 
        std::cout << "Please provide a positive integer to seed the PRNG." << std::endl;
        return 1;
    }

    int int_param = std::stoi(argv[1]); // give an integer to the CLI so that trivially parallelized nodes don't use same prng sequence
    int dim, S, init_L, seed, samples;
    bool random_start;
    RealType beta, t_max, log_base;

    std::ifstream infile("COPI_input.in");
    infile >> dim >> S >> init_L >> random_start >> beta >> seed >> samples >> t_max >> log_base;

    /*  // Pseudocode to run to convice yourself that the N-dimensional general screw-periodic boundary conditions are correct...
    for (int a = 0; a < sample.Total; ++a) {
        printf("neighb %d = ", a);
        for (int b = 0; b < 2 * sample.dim; ++b)
            printf("%d ", sample.Neighbors[a * (2 * sample.dim) + b]); 
        printf("\n");
    }  
    */
    
    // Route the runtime 'dim' variable to a compile-time template instantiation
    switch(dim) {
        case 2: // if dim is 2 (usually the case), call Run() with the template dim=2 and leave switch
            Run<2>(int_param, seed, init_L, beta, S, random_start, samples, t_max, log_base);
            break;
        case 3:
            Run<3>(int_param, seed, init_L, beta, S, random_start, samples, t_max, log_base);
            break;
        case 4:
            Run<4>(int_param, seed, init_L, beta, S, random_start, samples, t_max, log_base);
            break;
        default: // if dim isn't 2,3 or 4 in the input file, then run an error message
            std::cerr << "Error: Dimension " << dim << " is not supported. "
                      << "Please explicitly add 'case " << dim << ":' inside main()." << std::endl;
            return 1;
    }

    return 0;
}
