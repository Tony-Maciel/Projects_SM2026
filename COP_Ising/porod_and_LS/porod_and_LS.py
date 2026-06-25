import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.stats import linregress

# modify these 3 lines accordingly
L = 64
d = 4
file_pattern = "configs_64_0.600_*.dat"

N_sites = L**d
# Heuristics for the fitting window to extract the angular coefficient
# Skip early times to avoid initial quench transients 
t_min_fit = 1000000
# Stop fitting when the domain size hits the finite-size plateau
# Usually, finite size effects dominate when R(t) > L/2
R_max_fit = L / 2.0  

# SETUP RECIPROCAL SPACE (k-grid)
freqs = np.fft.fftfreq(L) * 2 * np.pi
freqs_shifted = np.fft.fftshift(freqs)
grids = np.meshgrid(*[freqs_shifted]*d, indexing='ij')
k_mag_grid = np.sqrt(sum(g**2 for g in grids))

dk = 2 * np.pi / L
max_k = np.pi
bins = np.arange(0, max_k + dk, dk)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Accumulators for the ensemble average
ensemble_S_k = {}   # dict mapping time -> sum of S(k) arrays
ensemble_counts = 0 # keep track of how many files we successfully read

# ENSEMBLE AVERAGING LOOP
file_list = glob.glob(file_pattern)
if not file_list:
    raise FileNotFoundError(f"No files matching '{file_pattern}' found in the directory.")

print(f"Found {len(file_list)} files. Beginning ensemble averaging...")

for file_idx, filename in enumerate(file_list):
    print(f"  -> Processing {filename} ({file_idx + 1}/{len(file_list)})")
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            data = np.fromstring(line, sep=' ')
            t = data[0]
            spin_flat = data[1:]
            
            if len(spin_flat) != N_sites:
                print(f"     [Warning] Expected {N_sites} spins at t={t}. Skipping.")
                continue
            
            spin_tensor = spin_flat.reshape([L]*d)
            
            # FFT
            fft_out = np.fft.fftn(spin_tensor)
            s_2d = (np.abs(fft_out)**2) / N_sites
            s_2d_shifted = np.fft.fftshift(s_2d)
            
            # Radial Binning
            hist_intensity, _ = np.histogram(k_mag_grid, bins=bins, weights=s_2d_shifted)
            hist_counts, _ = np.histogram(k_mag_grid, bins=bins)
            
            valid = hist_counts > 0
            k_pure = bin_centers[valid]
            s_k_pure = hist_intensity[valid] / hist_counts[valid]
            
            # Accumulate in the global dictionary
            if t not in ensemble_S_k:
                ensemble_S_k[t] = np.zeros_like(s_k_pure)
            ensemble_S_k[t] += s_k_pure
            
    ensemble_counts += 1

# CALCULATE OBSERVABLES FROM THE ENSEMBLE AVERAGE
print("\nCalculating characteristic length scales and extracting angular coefficient...")

times = np.array(sorted(ensemble_S_k.keys()))
R_t_list = []
S_k_avg_dict = {}

for t in times:
    # Average S(k) over all independent runs
    s_k_avg = ensemble_S_k[t] / ensemble_counts
    
    # Exclude the k=0 mode to avoid divergence and mean-magnetization artifacts
    k_pure = bin_centers[hist_counts > 0]
    mask = k_pure > 0 
    
    # Calculate the first moment k_1(t)
    k_1 = np.sum(k_pure[mask] * s_k_avg[mask]) / np.sum(s_k_avg[mask])
    R_t = 2 * np.pi / k_1
    
    R_t_list.append(R_t)
    S_k_avg_dict[t] = (k_pure, s_k_avg, k_1)

R_t_list = np.array(R_t_list)

# FIT THE ANGULAR COEFFICIENT (GROWTH EXPONENT)
# Filter data for the linear regression
fit_mask = (times >= t_min_fit) & (R_t_list <= R_max_fit)
fit_times = times[fit_mask]
fit_R_t = R_t_list[fit_mask]

if len(fit_times) > 1:
    # Perform linear regression on the log-log data
    log_t = np.log10(fit_times)
    log_R = np.log10(fit_R_t)
    slope, intercept, r_value, p_value, std_err = linregress(log_t, log_R)
    
    fit_line = (10**intercept) * (fit_times**slope)
    print(f"\n--- Fit Results ---")
    print(f"Fitted Angular Coefficient (Slope): {slope:.4f} ± {std_err:.4f}")
    print(f"Expected Lifshitz-Slyozov Value:  0.3333")
    print(f"R-squared:                        {r_value**2:.4f}")
    print(f"Time Window Used:                 {fit_times[0]} to {fit_times[-1]} MCS")
else:
    print("\nNot enough points in the valid scaling regime to perform a fit.")
    slope, intercept = 1/3, 0

# PLOTTING
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
fig.suptitle(f" {d}D COP Ising model with continuous time bulk diffusion algorithm ($L={L}$, {ensemble_counts} samples)", 
             fontsize=14, fontweight='bold')

# Lifshitz-Slyozov Scaling
ax1.loglog(times, R_t_list, 'ko-', label=r"$\langle R(t) \rangle$", markersize=5)

if len(fit_times) > 1:
    # Highlight the points used for the fit
    ax1.loglog(fit_times, fit_R_t, 'ro', markersize=5, 
               label="Data in Fit Window")
    # Plot the fitted line extending slightly past the data for visual clarity
    ext_t = np.array([fit_times[0]*0.001, fit_times[-1]*100.0])
    ext_R = (10**intercept) * (ext_t**slope)
    ax1.loglog(ext_t, ext_R, 'r--', linewidth=2, 
               label=f"Fit (Slope: {slope:.3f} +/- {std_err:.3f})")

# Plot theoretical expectation
t_ref = np.array([times[len(times)//4], times[-1]])
R_ref = R_t_list[len(times)//4] * (t_ref / t_ref[0])**(1/3)
ax1.loglog(t_ref, R_ref, 'b:', linewidth=2, label="Expected ($t^{1/3}$)")

ax1.set_xlabel("Time (a.u.)", fontsize=12)
ax1.set_ylabel(r"Characteristic domain size $R(t)$", fontsize=12)
# ax1.set_title("Dynamic Scaling & Angular Coefficient Fit", fontsize=12)
ax1.axhline(y=R_max_fit, color='gray', linestyle='-.', alpha=0.5, label=f"Plateau Cutoff ($L/2$)")
ax1.grid(True, which="both", ls="--", alpha=0.4)
ax1.legend()

# Dynamical Scaling Collapse 
plot_times = [times[i] for i in np.linspace(len(times)//4, len(times)-1, 6, dtype=int)]

for t in plot_times:
    k, s_k, k_1 = S_k_avg_dict[t]
    x = k / k_1
    F_x = s_k * (k_1**d)
    ax2.loglog(x, F_x, label=f"t = {int(t)}", alpha=0.8)

# Reference Porod tail for d=2 -> k^(-3)
porod_exponent = -(d + 1)
k_ref = np.array([2.0, 8.0])
S_ref = 3.0 * k_ref**(porod_exponent)
ax2.loglog(k_ref, S_ref, 'k--', linewidth=2.5, label=f"Porod's Law ($x^{{{porod_exponent}}}$)")

ax2.set_xlabel(r"Scaled wavevector $x = k / k_1(t)$", fontsize=12)
ax2.set_ylabel(r"$\langle S(k,t) \rangle \cdot k_1(t)^d$", fontsize=12)
ax2.set_title("Scaling Collapse & Porod Tail", fontsize=12)
ax2.set_xlim(0.1, 10) # Constrain view to ignore severe UV lattice artifacts
ax2.grid(True, which="both", ls="--", alpha=0.4)
ax2.legend()

plt.tight_layout()
plt.savefig(f"res_{L}_d{d}.png")
plt.show()
