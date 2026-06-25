In this subdirectory, I simulate brownian motion among an ensemble of particles in 2 dimensions starting
from the origin following the Wiener process: a spatial increment in a certain coordinate of a particle
is proportional to the square root of the time increment times a randomly generated variable from a normal
distribution with mean 0 and standard deviation 1.

The program brownian_motion.jl creates the output stored in the directories gifs/ and squared_displacement/. 

In gifs/ , there are some gifs showing the time evolution of a couple of particles and a circle whose radius 
grows proportionally to the square root of the current time of the simulation (this is supposed to guide the 
eyes as to the expected distance an average particle will be from the origin at any given time).

In squared_displacement/ , there are a couple plots which show the squared displacement averaged over an 
ensemble of particles over the course of time. The larger the ensemble, the more clearly a linear relationship 
develops, confirming theoretical predictions.

In dt_vs_sqrtdt/ , I test the statement proven in class that any other spatial increment that isn't proportional 
to the square root of the time increment yields incorrect dynamics. Namely, I test this with a spatial increment 
proportional to sqrt(dt) and just dt. As predicted in class, the increment with sqrt(dt) gives a realistic 
simulation of brownian motion (continuous analogue of the discrete random walker), while the increment with dt 
gives particles that have a squared displacement that grows linearly in time (so far so good...), BUT if dt is 
made smaller and smaller, the squared displacement grows slower in time (the rate at which the angular coefficient 
of the squared discplacement with respect to time decreases is linear), suggesting that as dt tends to zero, the 
particles get displaced less and less from the origin. In the limit of dt tending to zero, the particles simply
don't move from the origin, which corroborates the notion that a naive extrapolation to the continuous limit 
wouldn't give physical results.

The pdf Higham_numerical_SDE.pdf was used as supplementary material to making these simulations, hence why it's 
included.
