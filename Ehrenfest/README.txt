In this subdirectory, I simulate the Ehrenfest urn model up to a certain time t_max (e.g., 500) and 
average the number of balls in urn A over an ensemble of urns (e.g., 10, 100 etc.) at every given time.
Thus, at each time step, I will have an estimate of the mean number of balls in urn A (with error given 
by the standard deviation divided by sqrt(ensemble - 1)), which is used to compare against the analytical 
value of this quantity as derived in class.

For optimization purposes, I just store the urn A in memory and fix the total number of balls. So, 
instead of storing all the balls in a vector, randomly choosing one, seeing which urn it belongs to 
and swapping its urn, I just look at urn A and either remove a ball with probability N_A/N, or add a 
ball with complementary probability, where N is the total number of balls simulated (constant) and
N_A is the current number of balls in urn A. Clearly this yields the same dynamics as before.

The images in this directory are estimates of the mean number of balls in urn A compared with the 
expected analytical values for increasing values of samples in the ensemble. As expected, the larger 
the ensemble is, the closer the estimates come to the analytical values.
