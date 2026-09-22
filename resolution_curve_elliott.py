import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

E = np.array([5, 10, 20, 50, 100])
R = np.array([2.967, 2.028, 1.540, 0.995, 0.764]) / 100

def res(E, a, c):
    return np.sqrt(a**2/E + c**2)

popt, _ = curve_fit(res, E, R)
Efit = np.linspace(1, 105, 300)

plt.figure(figsize=(7, 5))
plt.plot(E, R, 'o', color='red')
plt.plot(Efit, res(Efit, *popt), color='blue')

plt.xlabel(r"$E$ [GeV]")
plt.ylabel(r"$\sigma_E/E$")
plt.xlim(0, 105)
plt.ylim(0.005, 0.07)

plt.savefig("plots/resolution_curve_elliott.png",
            dpi=300, bbox_inches="tight")

plt.show()

print(f"a = {popt[0]:.4f}")
print(f"c = {popt[1]:.4f}")
print("Saved: plots/resolution_curve_elliott.png")
