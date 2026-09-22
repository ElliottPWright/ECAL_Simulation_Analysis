import uproot
import awkward as ak
import numpy as np
import ROOT
from iminuit import Minuit
import glob
import os

# ============================================================
# INPUT / OUTPUT
# ============================================================

INPUT_PATH = "/eos/user/m/mahajanv/GeoModelSplitCal/elliott_sims/run_ECAL_Resolution/sim_mu-_*.root"

DATA_DIR = "data"
PLOT_DIR = "plots"

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)


# ============================================================
# LOAD MIP CALIBRATION
# ============================================================

calibration = np.load(os.path.join(DATA_DIR, "mip_calibration_elliott.npy"))

MIP_1CM = calibration[0]
MIP_6CM = calibration[1]

print("======================================")
print("MIP CALIBRATION")
print("======================================")
print(f"1 cm strips = {MIP_1CM:.6f} MeV/MIP")
print(f"6 cm strips = {MIP_6CM:.6f} MeV/MIP")
print()


# ============================================================
# GAUSSIAN FIT WITH IMINUIT
# ============================================================

def gaussian_fit(data, beam_energy):

    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    data = data[data > 0]

    if len(data) < 10:
        print("Not enough entries for Gaussian fit.")
        return None

    # Initial estimates
    mean0 = np.mean(data)
    sigma0 = np.std(data)

    # Restrict fit to a reasonable region around the response peak
    xmin = mean0 - 2.5 * sigma0
    xmax = mean0 + 2.5 * sigma0

    fit_data = data[(data >= xmin) & (data <= xmax)]

    if len(fit_data) < 10:
        print("Not enough entries inside Gaussian fit range.")
        return None

    # --------------------------------------------------------
    # Binned histogram
    # --------------------------------------------------------

    nbins = 100

    counts, edges = np.histogram(
        fit_data,
        bins=nbins,
        range=(xmin, xmax)
    )

    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    # Poisson uncertainties
    errors = np.sqrt(counts.astype(float))

    # Avoid zero-error bins
    errors[errors == 0] = 1.0

    # --------------------------------------------------------
    # Gaussian model
    # --------------------------------------------------------

    def gaussian(x, amplitude, mean, sigma):

        return (
            amplitude
            * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        )

    # Model amplitude is counts per bin
    amplitude0 = np.max(counts)

    # --------------------------------------------------------
    # Chi-square function
    # --------------------------------------------------------

    def chi2(amplitude, mean, sigma):

        expected = gaussian(
            centers,
            amplitude,
            mean,
            sigma
        )

        return np.sum(
            ((counts - expected) / errors) ** 2
        )

    # --------------------------------------------------------
    # iminuit
    # --------------------------------------------------------

    m = Minuit(
        chi2,
        amplitude=amplitude0,
        mean=mean0,
        sigma=sigma0
    )

    m.errordef = Minuit.LEAST_SQUARES

    m.limits["amplitude"] = (0, None)
    m.limits["sigma"] = (1e-6, None)
    m.limits["mean"] = (xmin, xmax)

    print("Running iminuit Gaussian minimization...")

    m.migrad()
    m.hesse()

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------

    fitted_mean = m.values["mean"]
    fitted_sigma = abs(m.values["sigma"])

    mean_error = m.errors["mean"]
    sigma_error = m.errors["sigma"]

    resolution = fitted_sigma / fitted_mean
    resolution_percent = 100.0 * resolution

    chi2_value = m.fval
    ndof = len(counts) - 3

    print()
    print(f"iminuit valid       = {m.valid}")
    print(f"chi2                = {chi2_value:.3f}")
    print(f"ndof                = {ndof}")
    print(f"chi2/ndof           = {chi2_value / ndof:.3f}")
    print(f"Mean                = {fitted_mean:.6f} MIP")
    print(f"Mean error          = {mean_error:.6f} MIP")
    print(f"Sigma               = {fitted_sigma:.6f} MIP")
    print(f"Sigma error         = {sigma_error:.6f} MIP")
    print(f"Resolution sigma/mu = {resolution:.6f}")
    print(f"Resolution          = {resolution_percent:.3f} %")

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    hist = ROOT.TH1D(
        f"h_{beam_energy}",
        f"{beam_energy} GeV electron response;Energy [MIP];Events",
        nbins,
        xmin,
        xmax
    )

    for value in fit_data:
        hist.Fill(float(value))

    canvas = ROOT.TCanvas(
        f"c_{beam_energy}",
        f"{beam_energy} GeV",
        900,
        700
    )

    hist.Draw("HIST")

    fit_function = ROOT.TF1(
        f"fit_{beam_energy}",
        "[0]*exp(-0.5*((x-[1])/[2])^2)",
        xmin,
        xmax
    )

    fit_function.SetParameters(
        m.values["amplitude"],
        fitted_mean,
        fitted_sigma
    )

    fit_function.SetLineWidth(2)
    fit_function.Draw("SAME")

    canvas.SaveAs(
        os.path.join(
            PLOT_DIR,
            f"electron_response_{beam_energy}GeV.png"
        )
    )

    return {
        "energy_GeV": beam_energy,
        "mean_MIP": fitted_mean,
        "mean_error_MIP": mean_error,
        "sigma_MIP": fitted_sigma,
        "sigma_error_MIP": sigma_error,
        "resolution": resolution,
        "resolution_percent": resolution_percent,
        "chi2": chi2_value,
        "ndof": ndof,
        "fit_valid": m.valid,
        "entries": len(data)
    }


# ============================================================
# ANALYZE ONE ELECTRON FILE
# ============================================================

def analyze_file(filename):

    print("======================================")
    print(f"FILE: {os.path.basename(filename)}")
    print("======================================")

    with uproot.open(filename) as f:

        tree = f["calo_events"]

        edep = tree["edep"].array(library="ak")
        typ = tree["type"].array(library="ak")

    event_energy_mip = []

    # --------------------------------------------------------
    # EVENTWISE ENERGY SUM
    # --------------------------------------------------------

    for event_edep, event_type in zip(edep, typ):

        energy_1cm = 0.0
        energy_6cm = 0.0

        for energy, strip_type in zip(
            event_edep,
            event_type
        ):

            energy = float(energy)
            strip_type = int(strip_type)

            # Thin PS strips = 1 cm
            if strip_type in (3, 4):
                energy_1cm += energy

            # Wide PVT strips = 6 cm
            elif strip_type in (1, 2):
                energy_6cm += energy

            # Types 5/6 = HPL/fibre
            # deliberately excluded

        # Convert each component to MIP units
        event_mip = (
            energy_1cm / MIP_1CM
            +
            energy_6cm / MIP_6CM
        )

        event_energy_mip.append(event_mip)

    event_energy_mip = np.asarray(
        event_energy_mip,
        dtype=float
    )

    # --------------------------------------------------------
    # Determine beam energy from filename
    # --------------------------------------------------------

    basename = os.path.basename(filename)

    beam_energy_mev = int(
        basename.split("_")[-2]
        .replace("MeV", "")
    )

    beam_energy_gev = beam_energy_mev / 1000.0

    print(f"Beam energy = {beam_energy_gev:g} GeV")
    print(f"Events      = {len(event_energy_mip)}")
    print(
        f"Mean raw    = {np.mean(event_energy_mip):.6f} MIP"
    )
    print(
        f"Std raw     = {np.std(event_energy_mip):.6f} MIP"
    )
    print()

    # Save eventwise distribution
    output_array = os.path.join(
        DATA_DIR,
        f"electron_{beam_energy_mev}MeV_MIP.npy"
    )

    np.save(
        output_array,
        event_energy_mip
    )

    # Gaussian fit
    result = gaussian_fit(
        event_energy_mip,
        beam_energy_gev
    )

    return result


# ============================================================
# MAIN
# ============================================================

files = sorted(glob.glob(INPUT_PATH))

print("Electron files found:")
for filename in files:
    print("  ", filename)

print()

results = []

for filename in files:

    try:

        result = analyze_file(filename)

        if result is not None:
            results.append(result)

    except Exception as e:

        print(
            f"ERROR processing {filename}: {e}"
        )

    print()


# ============================================================
# SAVE RESULTS
# ============================================================

if results:

    dtype = [
        ("energy_GeV", "f8"),
        ("mean_MIP", "f8"),
        ("mean_error_MIP", "f8"),
        ("sigma_MIP", "f8"),
        ("sigma_error_MIP", "f8"),
        ("resolution", "f8"),
        ("resolution_percent", "f8"),
        ("chi2", "f8"),
        ("ndof", "i4"),
        ("fit_valid", "?"),
        ("entries", "i4"),
    ]

    output = np.array(
        [
            (
                r["energy_GeV"],
                r["mean_MIP"],
                r["mean_error_MIP"],
                r["sigma_MIP"],
                r["sigma_error_MIP"],
                r["resolution"],
                r["resolution_percent"],
                r["chi2"],
                r["ndof"],
                r["fit_valid"],
                r["entries"],
            )
            for r in results
        ],
        dtype=dtype
    )

    output = np.sort(
        output,
        order="energy_GeV"
    )

    np.save(
        os.path.join(
            DATA_DIR,
            "resolution_results.npy"
        ),
        output
    )

    print("======================================")
    print("RESOLUTION RESULTS")
    print("======================================")

    for r in output:

        print(
            f"{r['energy_GeV']:6.0f} GeV : "
            f"mean = {r['mean_MIP']:.3f} MIP, "
            f"sigma = {r['sigma_MIP']:.3f} MIP, "
            f"resolution = {r['resolution_percent']:.3f} %"
        )

    print()
    print(
        "Saved: data/resolution_results.npy"
    )
