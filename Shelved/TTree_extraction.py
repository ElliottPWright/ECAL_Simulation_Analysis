 # ROOT File extraction executable

import os
import subprocess
#import ROOT
import csv
from scipy.signal import find_peaks, peak_widths
import shutil
import uproot
import scienceplots
import awkward as ak
import numpy as np                   
import matplotlib.pyplot as plt  
import argparse

#f = uproot.open(f'/eos/user/e/elwright/GeoModelSplitCal/build/calosim_out_10000MeV_mu-.root')                                                                                                                                      
plotting_directory = "/eos/user/e/elwright/EPWL_ECAL_Simulations/Figures"                                                                               
def extraction(file_input):                                        
    f= uproot.open(os.path.join("/eos/user/e/elwright/EPWL_ECAL_Simulations/output",file_input))
    tree = f["calo_events"]
    branches = ["edep", "x_global", "y_global"]

    data = tree.arrays(branches, library = "ak")

    Energies = ak.to_numpy(ak.flatten(data['edep']))

    counts, bin_edges = np.histogram(Energies, bins=50)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/ 2

    peak_index, _ = find_peaks(counts)
    max_peak = peak_index[np.argmax(counts[peak_index])]

    widths, width_heights, left_ips, right_ips = peak_widths(
            counts, [max_peak], real_height=0.5
            )

    fwhm = width[0] * (bin_centers[1] - bin_centers[0])

    #print(Energies)

    return Energies, fwhm, width_heights


def single_plotting(Energies, fwhm, width_heights, output):
    with plt.style.context(['science', 'no-latex']):
        plt.rcParams['figure.dpi'] = 200
        plt.hlines(*width_heights, color="C2")
        plt.hist(Energies, bins=50, histtype="step", color="blue")
        plt.xlabel("Energy [MeV]")
        plt.ylabel("Counts")
        plt.savefig((os.path.join(plotting_directory, output)))
        plt.show()

def plotting(fwhms, gaps, output):
        with plt.style.context(['science', 'no-latex']):
        plt.rcParams['figure.dpi'] = 200
        plt.plot(gaps, fwhm)
        plt.xlabel("Gap width [mm]")
        plt.ylabel("Resolution [MeV]")
        plt.savefig((os.path.join(plotting_directory, "SplitCalGapResolution.png")))
        plt.show()


def single_start():

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()

    file_input = args.input

    Energies, fwhm, width_heights = extraction(file_input)

    plotting(Energies, fwhm, width_heights, output=args.output)

def start():

    parser = argparse.ArgumentParser()

    Energies = []
    fwhms = []
    width_heights = []

    for file in os.listdir("/eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_test/"):
        Energy, fwhm, width_height = extraction(file)
        Energies.append(Energy)
        fwhms.append(fwhm)
        width_heights.append(width_height)


    plotting(fwhm, gaps, output=args.output)


if __name__ == "__main__":
    start()

