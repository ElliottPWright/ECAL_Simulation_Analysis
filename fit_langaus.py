#!/usr/bin/env python3
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def langaus_pdf(x, mpv, eta, sigma, A):
    std_x = (x - mpv) / eta
    landau = np.exp(-0.5 * (std_x + np.exp(-std_x))) / np.sqrt(2 * np.pi)
    gauss_term = np.exp(-0.5 * ((x - mpv) / sigma)**2)
    return A * (landau + 0.1 * gauss_term)

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

npy_files = sorted(glob.glob("data/edep_*MeV.npy"))

for npy_file in npy_files:
    energy_mev = int(re.search(r'(\d+)MeV', npy_file).group(1))
    energy_gev = energy_mev / 1000.0
    edep = np.load(npy_file)

    fig, ax = plt.subplots(figsize=(6, 4))
    counts, bin_edges, _ = ax.hist(edep, bins=100, histtype='step', color='black', label='Simulation Data')
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    peak_idx = np.argmax(counts)
    initial_mpv = bin_centers[peak_idx]
    
    try:
        popt_lg, _ = curve_fit(
            langaus_pdf, bin_centers, counts, 
            p0=[initial_mpv, 1.0, 2.0, np.max(counts)],
            maxfev=5000
        )
        x_fit = np.linspace(bin_centers[0], bin_centers[-1], 500)
        ax.plot(x_fit, langaus_pdf(x_fit, *popt_lg), 'r-', 
                label=f'LanGaus Fit\nMPV = {popt_lg[0]:.2f} MeV')
    except Exception as e:
        print(f"LanGaus fit failed for {energy_gev} GeV: {e}")

    ax.set_xlabel("E [MeV]")
    ax.set_ylabel("Counts")
    ax.set_yscale("log")
    ax.set_title(f"Raw Energy Deposition - {energy_gev:.0f} GeV")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"plots/langaus_fit_{energy_mev}MeV.png")
    plt.close()
    print(f"Saved: plots/langaus_fit_{energy_mev}MeV.png")
