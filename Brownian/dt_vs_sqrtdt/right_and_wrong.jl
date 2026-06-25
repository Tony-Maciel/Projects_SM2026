using Pkg 
Pkg.add("Plots")
Pkg.add("Random")
Pkg.add("Statistics")
using Random, Plots, Statistics

# Force the GR backend into headless mode (suppresses all visual output)
ENV["GKSwstype"] = "100"

const seed    = 1729    # To initialize PRNG (Xoshiro256++)
const t_max   = 100.0   # Total time simulated for a single system 
const samples = 1000     # Total number of samples 

"""
Simulates an ensemble of particles and returns the time array and the average squared displacement.
Accepts a boolean `use_correct_scaling` to toggle between dt and sqrt(dt).
"""
function simulate_ensemble(t_max, dt, samples, use_correct_scaling)
    Random.seed!(seed) 
    max_steps = round(Int, t_max / dt)
    
    # Toggle between the correct diffusion scaling and the incorrect linear scaling
    scale = use_correct_scaling ? sqrt(dt) : dt
    
    pos = zeros(Float64, samples, max_steps + 1, 2)
    
    for s in 1:samples
        for step in 2:(max_steps + 1)
            r1 = rand() 
            r2 = rand()
            # Box-Muller transformation
            pos[s, step, 1] = pos[s, step - 1, 1] + scale * sqrt(-2.0 * log(r1)) * cos(2pi * r2)
            pos[s, step, 2] = pos[s, step - 1, 2] + scale * sqrt(-2.0 * log(r1)) * sin(2pi * r2)
        end
    end
    
    time_array = 0:dt:t_max
    
    # Calculate r^2 = x^2 + y^2 for all samples and time steps
    r2 = pos[:, :, 1].^2 .+ pos[:, :, 2].^2 
    
    # Calculate the ensemble average and standard deviation over the first dimension (samples)
    avg_r2 = dropdims(mean(r2, dims=1), dims=1)
    std_r2 = dropdims(std(r2, dims=1), dims=1)
    
    # Standard error of the mean: std / sqrt(samples - 1)
    sem_r2 = std_r2 ./ sqrt(samples - 1)
    
    return time_array, avg_r2, sem_r2
end

function Run_Comparisons()
    dt_values = [0.2, 0.1, 0.05]
    
    # Plot incorrect dynamics (scale = dt)
    p_incorrect = plot(title="Incorrect Dynamics (Scale = dt)", 
                       xlabel="Time (t)", 
                       ylabel="Squared Displacement (r^2)", 
                       legend=:topleft,
                       linewidth=2)

    for dt in dt_values
        time_array, avg_r2, sem_r2 = simulate_ensemble(t_max, dt, samples, false)
        plot!(p_incorrect, time_array, avg_r2, ribbon=sem_r2, fillalpha=0.3, label="dt = $dt", linewidth=3)
    end
    
    savefig(p_incorrect, "incorrect_dynamics_comparison.png")

    # Plot correct dynamics (scale = sqrt(dt))
    p_correct = plot(title="Correct Dynamics (Scale = sqrt(dt))", 
                     xlabel="Time (t)", 
                     ylabel="Squared Displacement (r^2)", 
                     legend=:topleft,
                     linewidth=2)

    for dt in dt_values
        time_array, avg_r2, sem_r2 = simulate_ensemble(t_max, dt, samples, true)
        plot!(p_correct, time_array, avg_r2, ribbon=sem_r2, fillalpha=0.3, label="dt = $dt", linewidth=3, linestyle=:auto)
    end
    
    savefig(p_correct, "correct_dynamics_comparison.png")
    
    return nothing
end

Run_Comparisons()
