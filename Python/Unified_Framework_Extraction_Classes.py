# Code based on Matei's Event class in the ShipSoft framework
# booklet. Adjusted for the Mainz ECAL use case.

import ROOT
import os
import types
import sys
import numpy as np
import csv
from iminuit import Minuit
from scipy.optimize import curve_fit
import re
import glob
# import matplotlib.pyplot as plt

_alive = []

class SimFileReader:
    """
    This class opens and reads a ShipSoft simulation ROOT file.
    """

    def __init__(self, path, ntuple="events"):
        """
        Here we open a simulation file.
        """
    
        self.path = path
        self.reader = ROOT.RNTupleReader.Open(ntuple, path)
        self.views = {}
        _alive.append(self.reader)

    def view(self, field, cpp_type, optional=False):
            """Create a view on `field`; with optional=True return None if absent."""
            try:
                v = self.reader.GetView[cpp_type](field)
            except Exception:
                if optional:
                    return None
                raise
            _alive.append(v)
            self.views[field] = v
            return v
    
    def __getitem__(self, field):
        return self.views[field]
    
    def __len__(self):
        return int(self.reader.GetNEntries())
    
    def entries(self, limit=None):
        n = len(self)
        return range(n if limit is None else min(n, limit))



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
    
        Sim.view("event_number", "std::uint64_t")
        Sim.view("sim_hits", "std::vector<SHiP::SimHit>")
        Sim.view("sim_particles", "std::vector<SHiP::SimParticle>")

        self.z_positions = []
        self.edeps = []

        for i in Sim.entries():
            self.hits = Sim["sim_hits"](i)

            for h in self.hits:
                self.z_position = h.position[2]
                self.z_positions.append(self.z_position)

                if self.z_position < 500:
                    self.edep = h.energy_deposit
                    self.edeps.append(self.edep)

        self.edeps = 1000*np.array(self.edeps) # Conversion from GeV -> MeV

    def langaus_value(self, x: list, mpv: float, landau_sigma: float, gauss_sigma: float):
        """
        This method performs the Landau-Gaussian convolution.
        
        Parameters:
            x (array-like): x values.
            mpv (array-like): most-probable value, the distribution peak.
            landau_sigma, gauss_sigma (float): the distribution standard deviations.
        
        Returns:
            value: total * step. This represents one Langaus convolution value.
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

    def fit_mip(self, data: list, name: str, xmin: float, xmax: float):
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

        nbins = 10*int(xmax)

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


        def expected_counts(mpv: list, landau_sigma: float, gauss_sigma: float, amplitude: float):
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


        def nll(mpv: list, landau_sigma: float, gauss_sigma: float, amplitude: float):
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
        canvas.SetLeftMargin(0.15) 

        hist.GetXaxis().SetTitle("Energy deposited per strip [MeV]")
        hist.GetYaxis().SetTitle("Entries")

        # Change hist axis

        hist.GetXaxis().SetRangeUser(0, 30)
        hist.Draw()


        # ROOT function for visualisation only.
        # The minimisation itself was done by iminuit.

        def root_model(x, p):

            return (p[3] * self.langaus_value(x[0], p[0], p[1], p[2]) * bin_width)


        fit_func = ROOT.TF1("fit_" + name, root_model, xmin, xmax, 4)
        fit_func.SetParameters(mpv, landau_sigma, gauss_sigma, m.values["amplitude"])
        fit_func.SetLineColor(ROOT.kRed)
        fit_func.SetNpx(1000)
        fit_func.Draw("same")

        canvas.SaveAs(os.path.join(self.plotting_directory, rf"{name}.png"))

        return mpv, mpv_error

    def perform_fit(self, fit_min=1, fit_max=30):
        """
        Method to execute a MIP calibration on the strip layers.
        """

        mip, err = self.fit_mip(self.edeps, "MIP", fit_min, fit_max)
        
        output = np.array([mip, err])
        
        # Now save the values to a CSV
        
        with open(os.path.join(self.plotting_directory, 'MIP_Calibration.csv'), 'w', newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(["MIP", "Error"])
            writer.writerow(output)

        return output


class ECAL_Resolution(MIP_Calibration):
    """Determines the ECAL fractional resolution from sample electron data"""
    
    def __init__(self, input_directory: str, plotting_directory: str, calibration: list):
        """
        Initial method to collect the sum of energy deposits 
        across the events in a directory of simulation files.
        """
        self.plotting_directory = plotting_directory
    
        files = glob.glob(os.path.join(input_directory, "*.root"))
    
        def energy_key(name):
            m = re.search(r"(\d+(?:\.\d+)?)GeV", name)
            return float(m.group(1)) if m else float("inf")
    
        files.sort(key=energy_key)
    
        self.edep_sum_array = []

        for file in files:
            Sim = SimFileReader(file)
            Sim.view("event_number", "std::uint64_t")
            Sim.view("sim_hits", "std::vector<SHiP::SimHit>")
            Sim.view("sim_particles", "std::vector<SHiP::SimParticle>")
            
            self.z_positions = []
            self.edeps = []

            event_sums = []
            
            for i in Sim.entries():
                self.hits = Sim["sim_hits"](i)
                
                event_edeps = []
                
                for h in self.hits:
                    self.z_position = h.position[2]
                    self.z_positions.append(np.array(self.z_position))
            
                    if self.z_position < 500:
                        self.edep = h.energy_deposit
                        self.edeps.append(self.edep)
                        event_edeps.append(self.edep)
                
                
                edeps_MIP = [energy/calibration for energy in event_edeps]
                event_sums.append(np.sum(edeps_MIP))
                self.edep_sum_array.append(self.edep_sum)

            self.edep_sum_array.append(event_sums)

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


        def chi2(amplitude: float, mean: float, sigma: float):
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

        print("Running iminuit Gaussian minimisation...")

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

        canvas = ROOT.TCanvas(f"c_{beam_energy}", f"{beam_energy} GeV", 900, 700)

        hist.Draw("HIST")

        fit_function = ROOT.TF1(f"fit_{beam_energy}", "[0]*exp(-0.5*((x-[1])/[2])^2)", xmin, xmax)

        fit_function.SetParameters(m.values["amplitude"], fitted_mean, fitted_sigma)

        fit_function.SetLineWidth(2)
        fit_function.Draw("SAME")

        canvas.SaveAs(os.path.join(self.plotting_directory, f"electron_response_{beam_energy}GeV.png"))

        return {
            "energy_MeV": beam_energy,
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


    def perform_fit(self, beam_energy: list, plotting_directory: str):
        """
        Performs the energy resolution fit.

        Parameters:
            beam_energy (array-like): energy of the particles.
            plotting_directory (str): directory for plots.

        Returns:
            Plot of fractional energy resolution vs. beam energy
        """

        results = []

        for i in range(len(self.edep_sum_array)):
            result = self.gaussian_fit(self.edep_sum_array[i], beam_energy[i])
            results.append(result)

        results = np.array(results)

        dtype = [
            ("energy_MeV", "f8"),
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
                    r["energy_MeV"],
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

        columns = np.array(["energy_MeV", "mean_MIP", "mean_error_MIP", "resolution", "resolution_percent", "chi2", "ndof", "fit_valid", "entries"])

        with open(os.path.join(self.plotting_directory, 'ECAL_Resolution.csv'), 'w', newline='') as file:
                    writer = csv.writer(file, delimiter=',')
                    writer.writerow(columns)
                    for i in range(len(output)):
                        writer.writerow(output[i])


        def fractional_energy_resolution_fit(E, a, c):
            return np.sqrt(a**2/E + c**2)

        def resolution_fit_error(E: list, a: float, c: float, cov_matrix: tuple):
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


        popt, pcov = curve_fit(fractional_energy_resolution_fit, beam_energy, frac_resolutions)
        a, c = popt
        unc_params = np.sqrt(np.diag(pcov))

        fitted_energies = np.linspace(10**3, 10**5, num=10**3)
        fitted_frac_resolutions = fractional_energy_resolution_fit(fitted_energies, a, c)
        frac_error = resolution_fit_error(beam_energy, a, c, pcov)


        legend_text = (
            r"$\frac{\sigma_E}{E} = \sqrt{\frac{a^2}{E} + c^2}$" "\n"
            rf"$a = ({np.abs(a)/np.sqrt(1000)*100:.2f} \pm {unc_params[0]/np.sqrt(1000)*100:.2f}) \%$" + "\n"
            rf"$c = ({c*100:.2f} \pm {unc_params[1]*100:.2f}) \%$" "\n"
        )

        ## Need to remove the plt stuff here and replace it with ROOT canvas
        
        name = "ECAL_resolution"

        canvas = ROOT.TCanvas("c_" + name, name, 900, 700)

        graph = ROOT.TGraphErrors(len(beam_energy), beam_energy, frac_resolutions, 0*frac_error, frac_error)
        graph.GetXaxis().SetTitle("E [GeV]")
        graph.GetYaxis().SetTitle(r"$\frac{\sigma_E}{E}")
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(ROOT.kBlue)
        graph.SetName("SplitCalGap_energy_resolution")
        graph.Write()

        canvas.SaveAs(os.path.join(plotting_directory, "SplitCalGap_energy_resolution.png"))

        # with plt.style.context(['science', 'no-latex']):
        #     plt.rcParams['figure.dpi'] = 200
        #     plt.errorbar(beam_energy/10**3,frac_resolutions, yerr=frac_error, fmt='.', color='red', ecolor='black')
        #     plt.plot(fitted_energies/10**3, fitted_frac_resolutions, color='blue')
        #     plt.xlabel("E [GeV]")
        #     plt.ylabel(r"$\frac{\sigma_E}{E}$")
        #     plt.legend([legend_text], loc="best", prop={"family": "serif", "size": 11})
        #     plt.ylim([0, np.max(1.1*fitted_frac_resolutions)])
        #     plt.savefig(os.path.join(plotting_directory, "SplitCalGap_energy_resolution.png"))
        #     plt.show()

