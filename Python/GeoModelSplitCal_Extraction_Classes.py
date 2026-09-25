import os
import sys
sys.path.append(r"C:\root_v6.40.04\bin")
import uproot
import scienceplots
import matplotlib.pyplot as plt
import awkward as ak
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from iminuit import Minuit
from scipy.special import erf
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
import ROOT
import glob


class SimFileReader:
    """
    This class opens and reads a GeoModelSplitCal simulation ROOT file.
    """

    def __init__(self, file_path: str):
        """
        Here we open a simulation file.
        """

        self.file = uproot.open(file_path)

    def list_trees(self):
       """
       Here we return a list of all the trees available in the simulation file.
       """

       return [k for k in self.file.keys()]

    def get_branches(self, branches: list):
        """
        Get the relevant TTree branches from the ROOT file.
        """

        trees = self.list_trees()
        tree = self.file[trees[0]]  # Assuming the first tree is the one we want
        df = tree.arrays(branches, library="ak")
        data = [df[branch] for branch in branches]

        return data 

    def Diagnostic_Plotting(self, file_path: str, plotting_directory: str):
            """
            This method is for users to test if their simulations ran correctly.
            Currently, we plot: the full energy histogram; the first and last layer
            energy histogram; the full x-y hit distribution; the first and last layer
            hit distribution; and the energy-based acceptance.
            """

            print(f"Beginning diagnostics on: {file_path}")
            print(".............................................")

            # First of all, need to extract the relevant branches

            Sim = SimFileReader(os.path.join(file_path))
            branches = ["edep", "x_global", "y_global", "z_global", "layer"]
            edeps, x_globals, y_globals, z_globals, layers = Sim.get_branches(branches)

            # Now mask the arrays so that the relevant plots can be made

            first_layer_edeps = edeps[:][layers == 1]
            last_layer_edeps = edeps[:][layers == np.max(layers)]

            first_layer_edeps = ak.flatten(first_layer_edeps)
            last_layer_edeps = ak.flatten(last_layer_edeps)

            first_layer_x_hits = x_globals[:][layers == 1]
            last_layer_x_hits = x_globals[:][layers == np.max(layers)]
            first_layer_x_hits = ak.to_numpy(ak.flatten(first_layer_x_hits))
            last_layer_x_hits = ak.to_numpy(ak.flatten(last_layer_x_hits))
            
            first_layer_y_hits = y_globals[:][layers == 1]
            last_layer_y_hits = y_globals[:][layers == np.max(layers)]
            first_layer_y_hits = ak.to_numpy(ak.flatten(first_layer_y_hits))
            last_layer_y_hits = ak.to_numpy(ak.flatten(last_layer_y_hits))

            acceptance_energy = np.sum(last_layer_edeps)/np.sum(first_layer_edeps)

            # Checking that the sum of the energy histogram equals the input energy
            integrated_energy = np.array([np.sum(edeps[i]) for i in range(len(edeps))])
            integrated_energy_fraction = integrated_energy/10**4
            
            # print(f"Integrated energy = {integrated_energy:.2f}")
            # print(f"Fraction of original energy = {integrated_energy_fraction:.2f}")

            # print(f"Integrated energy = {integrated_energy}")
            # print(f"Fraction of original energy = {integrated_energy_fraction}")

            # All layers energy histogram
            print("Outputting all layers energy histogram.")

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(ak.to_numpy(ak.flatten(edeps)), bins = 100)
                plt.vlines(x=105.658, ymin=0, ymax=2*10**5, colors='red', linestyles='--', label='$\mu^{-}$ mass')
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.yscale("log")
                plt.legend()
                plt.savefig(os.path.join(plotting_directory, "SplitCal_all_layer_energy.png"))
                plt.show()
            

            # First layer energy histogram 
            print("Outputting first layer energy histogram")

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(first_layer_edeps, bins = 100)
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.yscale("log")
                plt.savefig(os.path.join(plotting_directory, "SplitCal_first_layer_energy.png"))
                plt.show()

            # Last layer energy histogram 
            print("Outputting last layer energy histogram")

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(last_layer_edeps, bins = 100)
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.yscale("log")
                plt.savefig(os.path.join(plotting_directory, "SplitCal_final_layer_energy.png"))
                plt.show()

            # x-z hits histogram
            print("Outputting all layers x-z hits histogram")
        
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(ak.to_numpy(ak.flatten(x_globals)), ak.to_numpy(ak.flatten(z_globals)), bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                plt.savefig(os.path.join(plotting_directory, "SplitCal_xz_hits.png"))
                plt.show()


            # y-z hits histogram
            print("Outputting all layers y-z hits histogram")
        
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(ak.to_numpy(ak.flatten(y_globals)), ak.to_numpy(ak.flatten(z_globals)), bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                plt.savefig(os.path.join(plotting_directory, "SplitCal_yz_hits.png"))
                plt.show()

            # All layers xy-hits histogram 
            print("Outputting all layers x-y hits histogram")
            
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(ak.to_numpy(ak.flatten(x_globals)), ak.to_numpy(ak.flatten(y_globals)), bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                plt.savefig(os.path.join(plotting_directory, "SplitCal_xy_all_hits.png"))
                plt.show()

            # First xy-hits histogram
            print("Outputting first layer x-y hits histogram")

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(first_layer_x_hits, first_layer_y_hits, bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                plt.savefig(os.path.join(plotting_directory, "SplitCal_xy_first_layer_hits.png"))
                plt.show() 

            # First layer xy-hits histogram
            print("Outputting last layer x-y hits histogram")

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(last_layer_x_hits, last_layer_y_hits, bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                plt.savefig(os.path.join(plotting_directory, "SplitCal_xy_last_layer_hits.png"))
                plt.show()
            
            # Acceptance energy-based
            print("Outputting energy-based acceptance")

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.plot(2160, acceptance_energy)
                plt.xlabel("Crack size [mm]")
                plt.ylabel("Acceptance")
                plt.savefig(os.path.join(plotting_directory, "SplitCal_energy_acceptance.png"))
                plt.show() 

    def close(self):
        """
        This simply closes the ROOT file you previously opened.
        """

        self.file.close()



class ProcessHits(SimFileReader):
    """
    Basic hits processing.
    """

    def __init__(self):
        self.first_layer_x_hits = self.x_globals[self.layers == 1]
        self.last_layer_x_hits = self.x_globals[self.layers == np.max(self.layers)]

        self.acceptance_energy = np.sum(self.last_layer_x_hits)/np.sum(self.first_layer_x_hits)
        self.acceptance_number = len(self.last_layer_x_hits)/len(self.first_layer_x_hits)



class MIP_Calibration(SimFileReader):
    """
    This class is dedicated to doing a minimum ionising particle
    (MIP) calibration.
    """

    def __init__(self, input_file: str, plotting_directory: str):
        """
        Here, we open the simulation file and extract the
        thin and wide layer energy depositions for later use
        """

        Sim = SimFileReader(input_file)

        self.plotting_directory = plotting_directory

        self.branches = ["edep", "type"]
        self.edeps, self.types = Sim.get_branches(self.branches)

        self.edeps = ak.flatten(self.edeps)
        self.types = ak.flatten(self.types)

        self.edep_wide = self.edeps[(self.types == 1) | (self.types == 2)]
        self.edep_thin = self.edeps[(self.types == 2) | (self.types == 3)]


    def langaus_value(self, x, mpv, landau_sigma, gauss_sigma):
        """
        This method performs the Landau-Gaussian convolution.

        Parameters:
            x (array-like): x values.
            mpv (array-like): most-probable value, the distribution peak.
            landau_sigma, gauss_sigma (float): the distribution standard deviations.
        
        Returns:
            value: total * step.
        """

        if landau_sigma <= 0 or gauss_sigma <= 0:
            return 0.0

        xmin = x - 5.0 * gauss_sigma
        xmax = x + 5.0 * gauss_sigma

        nsteps = 200
        step = (xmax - xmin) / nsteps

        total = 0.0

        for i in range(nsteps):

            xx_landau = xmin + (i + 0.5) * step

            landau = ROOT.TMath.Landau(xx_landau, mpv, landau_sigma, True)

            gaussian = ROOT.TMath.Gaus(x, xx_landau, gauss_sigma, True)

            total += landau * gaussian

        return total * step



    def fit_mip(self, data, name, xmin, xmax):
        """
        This method sets up the MIP calibration.

        Parameters:
            data (array-like): the simulation branch to be fitted.
            name (string): the name of the fit.
            xmin, xmax (floats): the maximum and minimum fittable domain.
        
        Returns:
            mpv, mpv_error (floats): distribution modal value and error.
        """


        # Start by creating the histogram 

        nbins = 100
        hist_min = 0.0
        hist_max = xmax

        hist = ROOT.TH1D(name, name, nbins, hist_min, hist_max)

        for value in data:
            hist.Fill(float(value))

        # Now fit bins strictly inside requested range

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


        # Finding the initial peak (mpv) based on the histogram

        peak_index = np.argmax(counts)
        peak = x_values[peak_index]
        max_content = counts[peak_index]

        print()
        print("======================================")
        print(name)
        print("======================================")
        print(f"Initial peak = {peak:.4f} MeV")
        print(f"Fit range    = {xmin:.2f} - {xmax:.2f} MeV")


        def expected_counts(mpv, landau_sigma, gauss_sigma, amplitude):
            """
            Here we find the expected counts based on the langaus model
            
            Parameters:
                mpv (array-like): most-probable value, the distribution peak.
                landau_sigma, gauss_sigma (float): the distribution standard deviations.
                amplitude (float): length of the data. 
            
            Returns:
                array: values

            """
            
            values = []

            for x in x_values:

                convolution = self.langaus_value(x, mpv, landau_sigma, gauss_sigma)

                # convolution is a probability density.
                # Multiply by bin width to obtain probability
                # contained in the histogram bin.
                expected = (amplitude * convolution * bin_width)

                values.append(expected)

            return np.asarray(values)


        def nll(mpv, landau_sigma, gauss_sigma, amplitude):
            """
            Extended Poisson negative log-likelihood.

            Parameters:
                mpv (array-like): most-probable value, the distribution peak.
                landau_sigma, gauss_sigma (float): the distribution standard deviations.
                amplitude (float): length of the data. 
            
            Returns:
                value: 2.0 * np.sum(term)
            """

            expected = expected_counts(mpv, landau_sigma, gauss_sigma, amplitude)

            # Protect against numerical problems
            expected = np.maximum(expected, 1e-12)

            # Poisson deviance
            term = np.where(counts > 0, expected - counts + counts * np.log(counts / expected), expected)

            return 2.0 * np.sum(term)


        # Running iminuit fit for landau-gaussian
        
        m = Minuit(nll, mpv=peak, landau_sigma=0.2, gauss_sigma=0.3, amplitude=len(data))

        m.limits["mpv"] = (1.8, 3.5)
        m.limits["landau_sigma"] = (0.05, 2.0)
        m.limits["gauss_sigma"] = (0.01, 1.0)
        m.limits["amplitude"] = (1.0, 1e7)

        # For -2 log likelihood
        m.errordef = Minuit.LIKELIHOOD

        print()
        print("Running iminuit minimisation...")

        m.migrad()

        if not m.valid:
            print("WARNING: MIGRAD did not converge.")
            print(m.fmin)

        m.hesse()

        # Printing out the results for the users to see

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

        # Creating and saving the relevant plots

        canvas = ROOT.TCanvas("c_" + name, name, 900, 700)

        hist.GetXaxis().SetTitle("Energy deposited per strip [MeV]")
        hist.GetYaxis().SetTitle("Entries")

        hist.Draw()

        # ROOT function for visualisation only.
        # The minimisation itself was done by iminuit.

        def root_model(x, p):

            return (p[3] * self.langaus_value(x[0], p[0], p[1], p[2]) * bin_width)

        fit_func = ROOT.TF1("fit_" + name, root_model, xmin, xmax, 4)

        fit_func.SetParameters(mpv, landau_sigma, gauss_sigma, m.values["amplitude"])

        fit_func.SetLineColor(ROOT.kRed)
        fit_func.Draw("same")

        canvas.SaveAs(os.path.join(self.plotting_directory, rf"{name}.png"))

        return mpv, mpv_error


    def perform_fit(self):
        """
        Method to execute a MIP calibration on the thin (1 cm)
        and wide (6 cm) strip layers.
        """
        
        mip_1cm, err_1cm = self.fit_mip(self.edep_thin, "MIP_1cm", np.min(self.edep_thin), np.max(self.edep_thin))

        mip_6cm, err_6cm = self.fit_mip(self.edep_wide, "MIP_6cm", np.min(self.edep_wide), np.max(self.edep_wide))

        return np.array([[mip_1cm, err_1cm], [mip_6cm, err_6cm]])



class ECAL_Resolution(MIP_Calibration):

    def __init__(self, input_directory: str, plotting_directory: str):
        """
        Initial method to collect the sum of energy deposits 
        across the events in a directory of simulation files.
        """
        self.plotting_directory = plotting_directory
        
        self.edep_sum_array = []

        for file in os.listdir(input_directory):
            Sim = SimFileReader(os.path.join(input_directory, file))
            self.branches = ["edep"]
            self.edeps = Sim.get_branches(self.branches)

            self.edep_sum = [np.sum(self.edeps[i]) for i in range(len(self.edeps))]
            self.edep_sum_array.append(self.edep_sum)

        self.edep_sum_array = np.array(self.edep_sum_array)


    def gaussian_fit(self, data: list, beam_energy: list):
        """
        Method to do a guassian fit with iminuit
        on the histogram of energy deposits across
        simulated events.

        Parameters:
            data (array-like): the simulation branch to be fitted.
            beam_energy (array-like): the beam energies used.
        
        Returns:
            data structure containing the fit parameters
        """
        
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

        # Producing the binned histogram

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


        def gaussian(x: list, amplitude: float, mean: float, sigma: float):
            """
            Gaussian fit.

            Parameters:
                x (array-like): x-data to be fitted.
                amplitude, mean, sigma (float): model parameters.

            Returns:
                Values of the guassian fit.
            """
            return (amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2))

        # Model amplitude is counts per bin
        amplitude0 = np.max(counts)


        def chi2(amplitude, mean, sigma):
            """
            Method to determine the chi-squared of a 
            fit.

            Parameters:
                amplitude, mean, sigma (float): model parameters.

            Returns:
                The chi-square value
            """
            
            expected = gaussian(centers, amplitude, mean, sigma)

            return np.sum(((counts - expected) / errors) ** 2)

        # iminuit fit

        m = Minuit(chi2, amplitude=amplitude0, mean=mean0, sigma=sigma0)
        m.errordef = Minuit.LEAST_SQUARES

        m.limits["amplitude"] = (0, None)
        m.limits["sigma"] = (1e-6, None)
        m.limits["mean"] = (xmin, xmax)

        print("Running iminuit Gaussian minimization...")

        m.migrad()
        m.hesse()


        # Printing the results of the iminuit fit for the user to see

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


        # Plotting the results of the iminuit for the user to see

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

        canvas.SaveAs(os.path.join(
                self.plotting_directory,
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

    
    def perform_fit(self, beam_energy):
        results = []
        
        for i in range(len(self.edep_sum_array)):
            result = self.gaussian_fit(self.edep_sum_array[i], beam_energy[i])
            results.append(result)

        results = np.array(results)

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

        output = np.array(output)

        def fractional_energy_resolution_fit(E, a, c):
            return np.sqrt(a**2/E + c**2)

        def resolution_fit_error(E, a, c, cov_matrix):
            """
            Calculate propagated errors for risetime derived from radii.

            Parameters:
                E (array-like): Measured E.
                a, c (float): Fit parameters.
                cov_matrix (2D array): Covariance matrix for the fit parameters [a, c].

            Returns:
                tuple: frac_res propagated uncertainties.
            """
            errors = []
            for Ei in E:
                d_y_da = a/(Ei*np.sqrt((a**2/Ei) + c**2))
                d_y_dc = c/np.sqrt((a**2/Ei) + c**2)

                # Gradient vector
                gradient = np.array([d_y_da, d_y_dc])
                
                # Propagated variance
                sigma_y_squared = gradient @ cov_matrix @ gradient.T
                sigma_y = np.sqrt(sigma_y_squared)
                errors.append(sigma_y)
                
            return np.array(errors)

        frac_resolutions = output["resolution"]


        popt, pcov = curve_fit(fractional_energy_resolution_fit, energies, frac_resolutions)
        a, c = popt
        unc_params = np.sqrt(np.diag(pcov))

        fitted_energies = np.linspace(10**3, 10**5, num=10**3)
        fitted_frac_resolutions = fractional_energy_resolution_fit(fitted_energies, a, c)
        frac_error = resolution_fit_error(energies, a, c, pcov)

        legend_text = (
            r"$\frac{\sigma_E}{E} = \sqrt{\frac{a^2}{E} + c^2}$" "\n"
            rf"$a = {a:.2f} \pm {unc_params[0]:.2f} \ [\sqrt{{\text{{GeV}}}}]$" + "\n"
            rf"$c = ({10**3*c:.2f} \pm {10**3*unc_params[1]:.2f})\times 10^{{-3}}$" "\n"
        )

        with plt.style.context(['science', 'no-latex']):
            plt.rcParams['figure.dpi'] = 200
            plt.errorbar(energies/10**3,frac_resolutions, yerr=frac_error, fmt='.', color='red', ecolor='black')
            plt.plot(fitted_energies/10**3, fitted_frac_resolutions, color='blue')
            plt.xlabel("E [GeV]")
            plt.ylabel(r"$\frac{\sigma_E}{E}$")
            plt.legend([legend_text], loc="best", prop={"family": "serif", "size": 11})
            #plt.yscale("log")
            #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
            plt.show()