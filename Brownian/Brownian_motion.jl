using Pkg 
Pkg.add("Plots")
Pkg.add("Random")
Pkg.add("Statistics")
# Pkg.add("BenchmarkTools")
using Random, Plots, Statistics
# using BenckmarkTools

# Force the GR backend into headless mode (suppresses all visual output)
ENV["GKSwstype"] = "100"

const seed      = 1729  # To initialize PRNG (Xoshiro256++)
const dt        = 0.01  # time discretization
const sqrtdt    = sqrt(dt)
const t_max     = 100.0 # Total time simulated a single system 
const samples   = 500   # Total number of samples simulated 
const max_steps = round(Int, t_max / dt)

"""Displays particles moving in time as gif and shows average squared distance over time."""
function Show(time_array, pos, t_max, samples)
    
    # Calculate r^2 = x^2 + y^2 for all samples and time steps
    r2 = pos[:, :, 1].^2 .+ pos[:, :, 2].^2 # .^2 is the broadcast operator in julia: it applies ^2 to all elements in an array

    # Calculate the ensemble average over the first dimension (samples)
    avg_r2 = dropdims(mean(r2, dims=1), dims=1)
    
    # Calculate the standard deviation across the samples
    std_r2 = dropdims(std(r2, dims=1), dims=1)
    
    # Calculate the standard error of the mean: std / sqrt(samples - 1)
    sem_r2 = std_r2 ./ sqrt(samples - 1)
    
    p_disp = plot(time_array, avg_r2, 
                  ribbon = sem_r2,           # Automatically creates the upper/lower bounds
                  fillalpha = 0.5,           # Transparency of the ribbon
                  fillcolor = :blue,         # Color of the ribbon
                  color = :red,              # Color of the average line
                  linewidth = 3, 
                  linestyle = :dot,
                  label = "Ensemble Average",
                  xlabel = "Time (t)", 
                  ylabel = "Squared Displacement (r^2)",
                  title = "Squared Displacement vs. Time ($samples samples)",
                  legend = :topleft)
                 
    # plot just a couple trajectories
    for s in 2:10
        plot!(p_disp, time_array, r2[s, :], label="", color=:gray, alpha=0.3)
    end
    
    savefig(p_disp, "squared_displacement_$samples.png")

    # gif...
    total_steps = length(time_array)
    frame_skip = max(1, div(total_steps, 200)) # Captures ~200 frames total
    
    # Determine static axis limits based on 3 standard deviations at t_max
    max_axis = 3 * sqrt(t_max)
    
    anim = @animate for step in 1:frame_skip:total_steps
        t = time_array[step]
        
        # Scatter plot of particle positions at this specific time step
        p_anim = scatter(pos[:, step, 1], pos[:, step, 2],
                         legend=false, markersize=3, markerstrokewidth=0,
                         color=:blue, alpha=0.7,
                         xlims=(-max_axis, max_axis), ylims=(-max_axis, max_axis),
                         aspect_ratio=:equal,
                         title="2D Brownian Motion (t = $(round(t, digits=1)))")
                         
        # The expanding variance circle (Radius = Standard Deviation = sqrt(t))
        # Using a 2-sigma boundary to encompass most particles (what really matters is the speed at which it moves)
        radius = 2*sqrt(t)
        theta = range(0, 2pi, length=100)
        circle_x = radius .* cos.(theta)
        circle_y = radius .* sin.(theta)
        
        plot!(p_anim, circle_x, circle_y, color=:red, linewidth=2, linestyle=:dash)
    end
    
    gif(anim, "brownian_particles_$samples.gif", fps=15)

    return nothing # prevents gif from being displayed
end

"""Simulates an ensemble of particles performing Brownian motion."""
function Run() 
    Random.seed!(seed) 
    
    # Initial condition: all particles start at origin at time zero
    pos = zeros(Float64, samples, max_steps + 1, 2) # 2 for the x and y dimensions
    
    for s in 1:samples
        for t in 2:(max_steps + 1) # time zero is index 1
            r1 = rand() 
            r2 = rand()
            # Box-Muller transformation
            # alternatively, could use randn() for built in function to generate from gaussian...
            pos[s, t, 1] = pos[s, t - 1, 1] + sqrtdt * sqrt(-2.0 * log(r1)) * cos(2pi * r2)
            pos[s, t, 2] = pos[s, t - 1, 2] + sqrtdt * sqrt(-2.0 * log(r1)) * sin(2pi * r2)
        end
    end
    
    time   = 0:dt:t_max
    Show(time, pos, t_max, samples)
    return nothing
end 

Run()

# replace above line with below if want to see how long it takes to run.
# (most of the time is spent including the libraries, hence why 10 and 100 samples takes roughly 3 seconds...)
# on my machine: 10 samples -> 3.215 s (1362062 allocations: 91.76 MiB)
#                100 samples -> 3.521 s (1524872 allocations: 131.60 MiB)
#                1000 samples -> 6.099 s (3350763 allocations: 524.80 MiB)
# @btime Run()
