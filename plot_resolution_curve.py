#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def resolution_func(E, a, b):
    return np.sqrt((a / np.sqrt(E))**2 + b**2)

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

data = np.load("data/fit_results.npy")
E_arr = data[:, 0]
mean_arr = data[:, 1]
sigma_arr = data[:, 2]

res_arr = sigma_arr / mean_arr

popt, pcov = curve_fit(resolution_func, E_arr, res_arr, p0=[0.10, 0.01])
a_fit, b_fit = popt

plt.figure(figsize=(6, 4))
plt.plot(E_arr, res_arr * 100, 'ko', label='Simulated Data')

E_smooth = np.linspace(min(E_arr) * 0.8, max(E_arr) * 1.1, 200)
plt.plot(E_smooth, resolution_func(E_smooth, a_fit, b_fit) * 100, 'r-', 
         label=f'Fit: $\sigma_E/E = {a_fit*100:.2f}\%/\sqrt{{E}} \oplus {b_fit*100:.2f}\%$')

plt.xlabel("Beam Energy $E$ [GeV]")
plt.ylabel("Energy Resolution $\sigma_E / E$ [%]")
plt.title("ECAL Energy Resolution Curve")
plt.grid(True, ls="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("plots/ecal_resolution_curve.png")

print(f"Stochastic term (a): {a_fit*100:.2f}%")
print(f"Constant term (b):   {b_fit*100:.2f}%")
