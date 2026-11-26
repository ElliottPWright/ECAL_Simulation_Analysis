from Unified_Framework_Extraction_Classes import SimFileReader, MIP_Calibration, ECAL_Resolution
import ROOT
import os
import types
import sys
import numpy as np

input_directory = r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/Unified_Framework/sim_runs_ECAL_Resolution/run_MIP_Calibration/sim_mu-_10GeV_1000mm.root"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/Unified_Framework/run_MIP_Calibration_3mm_Pb"


# MIP calibration

Sim = SimFileReader(input_directory)
MIP = MIP_Calibration(input_directory, plotting_directory)

calibration_3mm_Pb = MIP.perform_fit(fit_min=1, fit_max=10)

print()
print(".........................................")
print("Printing Results of MIP calibration")
print(calibration_3mm_Pb)


# # ECAL e- Resolution

energies = np.array([1, 3, 5, 7, 10, 15, 20, 30, 50, 100]) # Energies now in GeV

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/Unified_Framework/sim_runs_ECAL_Resolution/run_ECAL_Resolution"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/Unified_Framework/ECAL_Resolution_3mm_Pb"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_3mm_Pb[0])

results = ECAL.perform_fit(energies, plotting_directory)
