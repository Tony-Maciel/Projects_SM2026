using Pkg 
Pkg.add("Plots")
Pkg.add("Random")
Pkg.add("Statistics")
using Random, Plots, Statistics

const seed    = 1729  # To initialize PRNG (Xoshiro256++)
const N       = 100   # Total number of balls 
const t_max   = 500   # Total time simulated a single system 
const samples = 10000 # Total number of samples simulated 

"""Theoretical value for the average number of balls (total is N) in urn A at time t."""
function mu(N, t) 
  return N/2 + N/2 * (1.0 - 2.0/N)^t
end 

"""Displays average number of balls in urn A over the course of time."""
function Show(x, y, y_err, t_max, N)
    # Analytical values
    theoretical_y = [mu(N, t) for t in x]
    
    p = plot(x, y, ribbon=y_err, fillalpha=0.3, label="Simulated average", 
             xlabel="Time step (t)", ylabel="Balls in urn A",
             title="Ehrenfest Urn Model (N=100 and samples=10000)",
             linewidth=2, color=:blue)
             
    # Plot the theoretical expected values with a dotted line
    plot!(p, x, theoretical_y, label="Theoretical average", 
          linestyle=:dot, linewidth=3, color=:red)
          
    savefig(p, "mu$(samples)samples.png")
end 

"""Simulates an ensemble of Ehrenfest's urn model over a certain period of time."""
function Run() 
    Random.seed!(seed) 
    
    urn_A = zeros(Int, samples, t_max + 1)
    
    for s in 1:samples
        # Initial condition: all balls in urn A
        balls_in_A  = N
        urn_A[s, 1] = balls_in_A
        
        for t in 1:t_max
            # Probability of picking a ball from Urn A is (balls_in_A / N)
            if rand() < (balls_in_A / N)
                balls_in_A -= 1 # Ball moves from A to B
            else
                balls_in_A += 1 # Ball moves from B to A
            end
            urn_A[s, t + 1] = balls_in_A
        end
    end
    
    # Average over all ensembles at each time step
    time   = 0:t_max
    avg_mu = vec(mean(urn_A, dims=1)) # dims=1 means to average over the columns
    err_mu = vec(std(urn_A, dims=1) / sqrt(samples - 1)) 
    
    Show(time, avg_mu, err_mu, t_max, N)
end 

Run()
