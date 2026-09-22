#!/usr/bin/env python3

import os
import numpy as np
import ROOT
from iminuit import Minuit

os.makedirs("plots", exist_ok=True)
os.makedirs("data", exist_ok=True)

# ============================================================
# Load muon strip-energy data
# ============================================================

# ============================================================
# Load Elliott's muon ROOT files
# ============================================================

import uproot
import awkward as ak

MIP_FILES = [
    "/eos/user/m/mahajanv/GeoModelSplitCal/elliott_sims/run_MIP_Calibration/sim_mu-_10000MeV_1000mm.root",
    "/eos/user/m/mahajanv/GeoModelSplitCal/elliott_sims/run_MIP_Calibration/sim_mu-_50000MeV_1000mm.root",
]

data_1cm_list = []
data_6cm_list = []

for filename in MIP_FILES:

    print(f"Reading: {filename}")

    tree = uproot.open(filename)["calo_events"]

    edep = tree["edep"].array(library="ak")
    types = tree["type"].array(library="ak")

    # Loop over events
    for event_edep, event_type in zip(edep, types):

        # 1 cm strips = types 3, 4
        mask_1cm = (event_type == 3) | (event_type == 4)

        # 6 cm strips = types 1, 2
        mask_6cm = (event_type == 1) | (event_type == 2)

        data_1cm_list.extend(
            np.asarray(event_edep[mask_1cm], dtype=float)
        )

        data_6cm_list.extend(
            np.asarray(event_edep[mask_6cm], dtype=float)
        )

data_1cm = np.asarray(data_1cm_list)
data_6cm = np.asarray(data_6cm_list)

# Keep only non-zero deposits
data_1cm = data_1cm[data_1cm > 0]
data_6cm = data_6cm[data_6cm > 0]

print(f"1 cm entries: {len(data_1cm)}")
print(f"6 cm entries: {len(data_6cm)}")

# ============================================================
# Landau-Gaussian convolution
# ============================================================

def langaus_value(x, mpv, landau_sigma, gauss_sigma):

    if landau_sigma <= 0 or gauss_sigma <= 0:
        return 0.0

    xmin = x - 5.0 * gauss_sigma
    xmax = x + 5.0 * gauss_sigma

    nsteps = 200
    step = (xmax - xmin) / nsteps

    total = 0.0

    for i in range(nsteps):

        xx_landau = xmin + (i + 0.5) * step

        landau = ROOT.TMath.Landau(
            xx_landau,
            mpv,
            landau_sigma,
            True
        )

        gaussian = ROOT.TMath.Gaus(
            x,
            xx_landau,
            gauss_sigma,
            True
        )

        total += landau * gaussian

    return total * step


# ============================================================
# Fit
# ============================================================

def fit_mip(data, name, xmin, xmax):

    nbins = 100
    hist_min = 0.0
    hist_max = xmax

    hist = ROOT.TH1D(
        name,
        name,
        nbins,
        hist_min,
        hist_max
    )

    for value in data:
        hist.Fill(float(value))

    # --------------------------------------------------------
    # Fit bins strictly inside requested range
    # --------------------------------------------------------

    x_values = []
    counts = []

    bin_width = hist.GetBinWidth(1)

    for i in range(1, nbins + 1):

        x = hist.GetBinCenter(i)

        if x < xmin or x > xmax:
            continue

        y = hist.GetBinContent(i)

        if y <= 0:
            continue

        x_values.append(x)
        counts.append(y)

    x_values = np.asarray(x_values)
    counts = np.asarray(counts)

    # --------------------------------------------------------
    # Initial peak
    # --------------------------------------------------------

    peak_index = np.argmax(counts)
    peak = x_values[peak_index]
    max_content = counts[peak_index]

    print()
    print("======================================")
    print(name)
    print("======================================")
    print(f"Initial peak = {peak:.4f} MeV")
    print(f"Fit range    = {xmin:.2f} - {xmax:.2f} MeV")

    # --------------------------------------------------------
    # Model: expected histogram counts
    # --------------------------------------------------------

    def expected_counts(
        mpv,
        landau_sigma,
        gauss_sigma,
        amplitude
    ):

        values = []

        for x in x_values:

            convolution = langaus_value(
                x,
                mpv,
                landau_sigma,
                gauss_sigma
            )

            # convolution is a probability density.
            # Multiply by bin width to obtain probability
            # contained in the histogram bin.
            expected = (
                amplitude
                * convolution
                * bin_width
            )

            values.append(expected)

        return np.asarray(values)

    # --------------------------------------------------------
    # Extended Poisson negative log-likelihood
    # --------------------------------------------------------

    def nll(
        mpv,
        landau_sigma,
        gauss_sigma,
        amplitude
    ):

        expected = expected_counts(
            mpv,
            landau_sigma,
            gauss_sigma,
            amplitude
        )

        # Protect against numerical problems
        expected = np.maximum(expected, 1e-12)

        # Poisson deviance
        term = np.where(
            counts > 0,
            expected - counts + counts * np.log(counts / expected),
            expected
        )

        return 2.0 * np.sum(term)

    # --------------------------------------------------------
    # iminuit
    # --------------------------------------------------------

    m = Minuit(
        nll,
        mpv=peak,
        landau_sigma=0.2,
        gauss_sigma=0.3,
        amplitude=len(data)
    )

    m.limits["mpv"] = (1.8, 3.5)
    m.limits["landau_sigma"] = (0.05, 2.0)
    m.limits["gauss_sigma"] = (0.01, 1.0)
    m.limits["amplitude"] = (1.0, 1e7)

    # For -2 log likelihood
    m.errordef = Minuit.LIKELIHOOD

    print()
    print("Running iminuit minimization...")

    m.migrad()

    if not m.valid:
        print("WARNING: MIGRAD did not converge.")
        print(m.fmin)

    m.hesse()

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------

    mpv = m.values["mpv"]
    mpv_error = m.errors["mpv"]

    landau_sigma = m.values["landau_sigma"]
    gauss_sigma = m.values["gauss_sigma"]

    print()
    print(f"iminuit valid = {m.valid}")
    print(f"NLL/deviance   = {m.fval:.3f}")

    print()
    print(f"MPV            = {mpv:.6f} MeV")
    print(f"MPV error      = {mpv_error:.6f} MeV")
    print(f"Landau sigma   = {landau_sigma:.6f} MeV")
    print(f"Gaussian sigma = {gauss_sigma:.6f} MeV")
    print(f"Amplitude      = {m.values['amplitude']:.3f}")

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    canvas = ROOT.TCanvas(
        "c_" + name,
        name,
        900,
        700
    )

    hist.GetXaxis().SetTitle(
        "Energy deposited per strip [MeV]"
    )

    hist.GetYaxis().SetTitle(
        "Entries"
    )

    hist.Draw()

    # ROOT function for visualisation only.
    # The minimisation itself was done by iminuit.

    def root_model(x, p):

        return (
            p[3]
            * langaus_value(
                x[0],
                p[0],
                p[1],
                p[2]
            )
            * bin_width
        )

    fit_func = ROOT.TF1(
        "fit_" + name,
        root_model,
        xmin,
        xmax,
        4
    )

    fit_func.SetParameters(
        mpv,
        landau_sigma,
        gauss_sigma,
        m.values["amplitude"]
    )

    fit_func.SetLineColor(ROOT.kRed)
    fit_func.Draw("same")

    canvas.SaveAs(
        f"plots/{name}.png"
    )

    return mpv, mpv_error


# ============================================================
# Perform fits
# ============================================================

mip_1cm, err_1cm = fit_mip(
    data_1cm,
    "MIP_1cm",
    1.8,
    3.5
)

mip_6cm, err_6cm = fit_mip(
    data_6cm,
    "MIP_6cm",
    1.8,
    3.5
)


# ============================================================
# Save calibration
# ============================================================

calibration = np.array([
    mip_1cm,
    mip_6cm
])

np.save(
    "data/mip_calibration_elliott.npy",
    calibration
)

print()
print("======================================")
print("MIP CALIBRATION")
print("======================================")

print(
    f"1 cm strips = "
    f"{mip_1cm:.6f} MeV/MIP"
)

print(
    f"6 cm strips = "
    f"{mip_6cm:.6f} MeV/MIP"
)

print()
print("Saved:")
print("data/mip_calibration.npy")
