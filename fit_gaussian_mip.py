#!/usr/bin/env python3
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-0.5 * ((x - mean) / sigma)**2)

# Built-in clean styling
plt.style.use('seaborn-v0_8-ticks')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'figure.dpi': 200
})

MIP_SCALE_FACTOR = 0.211  # MeV per MIP
npy_files = sorted(glob.glob("data/edep_*MeV.npy"))
fit_results = []

for npy_file in npy_files:
    energy_mev = int(re.search(r'(\d+)MeV', npy_file).group(1))
    energy_gev = energy_mev / 1000.0
    edep = np.load(npy_file)
    
    edep_mip = edep / MIP_SCALE_FACTOR

    fig, ax = plt.subplots(figsize=(6, 4))
    counts, bin_edges, _ = ax.hist(edep_mip, bins=100, histtype='step', color='black', label='MIP Calibrated')
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    mean_init = np.mean(edep_mip)
    sigma_init = np.std(edep_mip)

    try:
        popt_g, _ = curve_fit(
            gaussian, bin_centers, counts, 
            p0=[np.max(counts), mean_init, sigma_init]
        )
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 500)
        
        mean_fit = popt_g[1]
        sigma_fit = abs(popt_g[2])
        
        ax.plot(x_fit, gaussian(x_fit, *popt_g), 'r--', 
                label=f'Gaussian Fit\n$\mu$ = {mean_fit:.1f} MIPs\n$\sigma$ = {sigma_fit:.1f} MIPs')
        
        fit_results.append((energy_gev, mean_fit, sigma_fit))
    except Exception as e:
        print(f"Gaussian fit failed for {energy_gev} GeV: {e}")

    ax.set_xlabel("E [MIP]")
    ax.set_ylabel("Counts")
    ax.set_yscale("log")
    ax.set_title(f"MIP Calibrated Response - {energy_gev:.0f} GeV")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"plots/gaussian_fit_mip_{energy_mev}MeV.png")
    plt.close()

np.save("data/fit_results.npy", np.array(fit_results))
print("Saved fit results -> data/fit_results.npy")
