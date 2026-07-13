# This project contains many programs that simulate some models studied in my statistical mechanics classes.
  ## 1) Ehrenfest urn model (entree). 
  The Ehrenfest urn model can be stated as follows: there are two urns, A and B, that both initially contain a certain number of identical balls.
  In total, there are N (identical) balls. At each time step, a ball must randomly be chosen, and the urn with which it belongs to must be swapped.
  After a sufficiently long time, it's expected that both urns will have the same number of balls on average.

  To illustrate, suppose all the N balls start in urn A. Then in the first time step, a ball is randomly chosen. Since all the balls are in urn A, 
  this chosen ball will belong to urn A. Then this ball will get moved to urn B. In the next time step, another ball is randomly chosen, but this 
  time a ball from urn B could be chosen (with probability 1/N), and this whole process is repeated ad nauseum.
  
  ## 2) Conserved order parameter Ising model with the enhanced bulk diffusion algorithm (main course).
  The conserved order parameter Ising model ("COP Ising model") is just the ordinary Ising model (a cubic lattice with sites with 2 states in each
  vertex that can only interact with its nearest neighbors, and has an interaction energy associated with each nearest neighbor bond), but the 
  total magnetization (order parameter), is fixed at a certain value. 

  Because of this restriction on the magnetization, conventional monte carlo algorithms won't work. For example, in the Metropolis algorithm, a 
  spin could be flipped (or not), which results in a different magnetization (not allowed). To get around this, many algorithms have been proposed
  that satisfy ergodicity and detailed balance just like the Metropolis algorithm (so the correct Boltzmann distribution is sampled), but don't
  alter the magnetization. This can be done by, for example, swapping the states of 2 spins, instead of just flipping one.

  If a low temperature is simulated for this system, it's expected that large magnetic domains will form, since the free energy is minimized by
  minimizing the surface area of these magnetic domains, i.e., gathering all the smaller domains into one big one.

  The most common algorithm for this model is the Kawasaki algorithm (imagine the Metropolis algorithm, but instead of flipping a spin, it 
  swaps a spin with its nearest neighbor). But, since I want to learn an algorithm that is as different as possible from the ones I already
  know, I decided to use the enhanced bulk diffusion algorithm (c.f., Newman and Barkema), which also satisfies ergodicity and detailed 
  balance, but gives a higher probability for spins to be swapped among domains (when compared with the Kawasaki algorithm), facilitating 
  the formation of large magnetic domains.
  
  ## 3) Brownian motion (dessert). 
  Brownian motion is the random motion of particles inside some fluid. Mathematically, it's described as a Wiener process, which can be 
  stated as a couple of axioms, but, in my opinion, isn't obvious at all that it describes the erratic little dance that particles do 
  when suspended in water, for example. 

  I believe that these axioms were probably motivated by analyzing the divergent limits that appear when trying to turn the discrete random 
  walker problem into a continuous model (both in space and time). When trying to do this, because of the central limit theorem, you arrive 
  at the natural conclusion that after each (infinitesimal) time step, the (infinitesimal) spatial step must be incremented by an amount 
  proportional to the SQUARE ROOT of the time step, NOT by an amount proportional to the time step. If this wasn't the case, then there 
  wouldn't be a diffusion process at all! I verify this by trying both approaches and seeing that only with the square root increment of 
  the time step, the motion of the particles actually resemble that of ordinary diffusive motion (but in the continuous limit).
