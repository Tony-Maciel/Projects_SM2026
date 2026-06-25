In this directory, COPI.cpp simulates the conserved order parameter Ising model in 2 dimensions with 
periodic boundary conditions. However, it's programmed to be able to handle generalized helical boundary 
conditions and arbitrary dimensions (as long as the latter is included in the end of main() for optimization
reasons). The algorithm used is the same one proposed by Marko and Barkema (c.f. COP_IsingModel-PhysRevE.52.2522.pdf).
This algorithm is also explained in p. 114 of S. Puri and V. Wadhawan (Kinetics_of_phase_transitions.pdf).
It's a continuous time, rejection-free algorithm that is very similar to the Gillespie algorithm, hence
why it was used. However, it's known that it isn't the *most* efficient algorithm in simulating this model
for some values of temperatures (like the one I decided to simulate... :|), as pointed out by Newman and Barkema, p.282.

(Bray_theory_of_phase_ordering_kinetics.pdf contains the dynamical scaling hypothesis used in corroborating Porod's law 
(c.f., eq. 7).)

This directory uses C++ program to do the actual simulations themselves and python files for the data analysis and 
visualization. Standard practice dictates that a virtual environment should be created in order to use python files. 
So if you want to run any python program for yourself, run: 

python -m venv myenv source myenv/bin/activate pip install 

This creates a virtual environment. After this, anytime you want to run a python program in this directory, you should run: 

source myenv/bin/activate 

Then try running the python files. At first they will give errors saying that you don't have certain libraries, so you should 
simply pip install to resolve any dependencies. Afterwards, if you source this directory again, you won't have to pip install 
from scratch. All the libraries used were: numpy, matplotlib, glob and scipy.

The gist of the C++ program is as follows: it reads COPI_input.in to get the parameters used for the simulation (e.g., dimension, 
system size L, inverse temperature beta etc.), then it runs multiple simulations in parallel until it has simulated the 
desired number of systems, all the while saving to a unique output file the current simulation time and the configuration of 
the system at that instant (this is done in logarithmic spacing to avoid huge output files, and because of convenience for plotting
later). After this program finishes running, the desired python file should have its parameters modified to read the desired 
output file(s) and generate the correct results (inside the respective directory).

Inside the directory gifs/ , there are some gifs showing snapshots of this model for a couple system sizes as well as the 
python file and output files that were used to generate these gifs, which are inside a .zip file.

All the data in the directories here are in .zip files because they are pretty large (the total size of all data files is
around 111Gb (:O). But around 100Gb of this is just due to the simulations on the 4D lattice with 64^4 spins, so unpack
that zip file with caution!)

Inside the directory porod_and_LS/ , there are images corroborating the: Lifshitz-Slyozov law, which states that the average 
size of a magnetic domain grows like t^(1/3) in the conserved order parameter Ising model (c.f. Newman and Barkema, p. 269);
and Porod's law, which essentially states that, after a given amount of time, sharp interfaces will form, corresponding to 
the magnetic domains which grow larger and larger with time. Mathematically, Porod's law states that as the magnitude of the 
wavevector k tends to infinity, the circularly averaged structure factor tends to zero as k^(-(dim + 1)). 

As for some of the optimizations used in the C++ file, here I will name a couple of the more important ones:
  1)A template is used for certain values of dim (dim=2,3 and 4. Any other dimension apart from these is uncommon) to allow
    the compiler to aggressively optimize further. This entails loop unrolling, statically allocated arrays etc. Because of 
    this optimization, I was able to simulate systems with roughly 17 million spins in 4 dimensions in roughly the same amount 
    of time as it took to simulate systems with roughly 2 million spins in 3 dimensions! This is because a lot of data 
    structures that were constantly used had sizes proportional to the dimension of the system. 

  2)openMP is used to simulate the samples in parallel, which gives a roughly linear increase in speedup (since this task is 
    trivially parallelizable). To facilitate this parallel task, I created structs (System) and instantiated them for each 
    thread.

  3)Multiple data structures were created to optimize the simulation, such as: a vector that precomputes exponentials before 
    the simulation is actually carried out (this is only possible because each spin can assume discrete values), a vector 
    that contains the addresses of each site's neighbors (this avoids calculating modulo function calls ad hoc, and is 
    different depending on the boundary conditions desired), a vector that contains the addresses of each neighbor that should 
    have its information updated after a pair of spins has its spins flipped etc.

  4)Xoshiro256+ was used as the pseudo random number generator. It is one of the fastest prngs today, generating a single 
    random number in under a nanosecond for most modern computers, and it passes all known statistical tests (according to 
    the authors of this prng). For every serious monte carlo simulation, the prng used is almost always the bottleneck of 
    the code, so even just a 5% speedup in a prng can generate considerably faster code. (alternatively, if you want the 
    highest quality random numbers possible, then the gold standard is Ranlux and Ranlux++. In my experience, a monte carlo 
    program that uses these prngs will take longer to run, and in the end it won't be worth it, because in that same amount 
    of time, I could have run my simulations with Xoshiro256+, simulating larger lattices or more temperature or more monte 
    carlo steps, and gotten better results. This isn't to say that Ranlux doesn't have its place of course, because a lot 
    of scientists, especially those in high energy physics, scrupulously adhere to using only Ranlux, because according to 
    them, the quality of the prng really does make a difference.)
