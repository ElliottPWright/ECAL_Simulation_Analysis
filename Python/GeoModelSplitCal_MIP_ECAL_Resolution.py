## Imports

import os
import matplotlib.pyplot as plt
import awkward as ak
import scienceplots
import numpy as np
from GeoModelSplitCal_Extraction_Classes import SimFileReader, ProcessHits, MIP_Calibration, ECAL_Resolution


## 3 mm Pb Hextant Centred

# MIP calibration

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_MIP_Calibration_Beam_Centred"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/run_MIP_Calibration_3mm_Pb_Centred_Beam"

Sim = SimFileReader(os.path.join(file_path, "sim_mu-_10000MeV_1000mm.root"))
MIP = MIP_Calibration(os.path.join(file_path, "sim_mu-_10000MeV_1000mm.root"), plotting_directory)

calibration_3mm_Pb_Centred = MIP.perform_fit()

# ECAL Resolution

energies = np.array([1_000, 3_000, 5_000, 7_000, 10_000, 15_000, 20_000, 30_000, 50_000, 100_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_ECAL_Resolution_Beam_Centred"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_3mm_Pb_Centred_Beam"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_3mm_Pb_Centred[:,0])
ECAL.perform_fit(energies, plotting_directory)


## 3 mm Pb Strip Centred

# MIP calibration

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_MIP_Calibration"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/run_MIP_Calibration_3mm_Pb"

Sim = SimFileReader(os.path.join(file_path, "sim_mu-_10000MeV_1000mm.root"))
MIP = MIP_Calibration(os.path.join(file_path, "sim_mu-_10000MeV_1000mm.root"), plotting_directory)

calibration_3mm_Pb = MIP.perform_fit()

# ECAL e- Resolution

energies = np.array([1_000, 3_000, 5_000, 7_000, 10_000, 15_000, 20_000, 30_000, 50_000, 100_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_ECAL_Resolution"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_3mm_Pb"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_3mm_Pb[:,0])
ECAL.perform_fit(energies, plotting_directory)

# ECAL gamma Resolution

energies = np.array([1_000, 3_000, 5_000, 7_000, 10_000, 15_000, 20_000, 30_000, 50_000, 100_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_ECAL_Resolution_Gamma"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_Gamma_3mm_Pb"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_3mm_Pb[:,0])
ECAL.perform_fit(energies, plotting_directory)



## 9 mm Fe Strip Centred

# MIP calibration (Vani's)

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/vmahajan/fe_mu_runs"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/run_MIP_Calibration_9mm_Fe_VANI"

Sim = SimFileReader(os.path.join(file_path, "sim_mu-_10k_10000MeV.root"))
MIP = MIP_Calibration(os.path.join(file_path, "sim_mu-_10k_10000MeV.root"), plotting_directory)


calibration_9mm_Fe_VANI = MIP.perform_fit()

# ECAL Resolution (Vani's)

energies = np.array([5_000, 10_000, 20_000, 50_000, 100_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/vmahajan/fe_e-_runs"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_9mm_Fe_VANI"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_9mm_Fe_VANI[:,0])
ECAL.perform_fit(energies, plotting_directory)


# MIP calibration (Elliott's)

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_MIP_Calibration_9mm_Fe"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/run_MIP_Calibration_9mm_Fe"

Sim = SimFileReader(os.path.join(file_path, "sim_mu-_10k_10000MeV.root"))
MIP = MIP_Calibration(os.path.join(file_path, "sim_mu-_10k_10000MeV.root"), plotting_directory)


calibration_9mm_Fe = MIP.perform_fit()



# ECAL e- Resolution (Elliott's)

energies = np.array([5_000, 10_000, 20_000, 50_000, 100_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_ECAL_Resolution_9mm_Fe"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_9mm_Fe"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_9mm_Fe[:,0])
ECAL.perform_fit(energies, plotting_directory)


# ECAL gamma Resolution (Elliott's)

energies = np.array([5_000, 10_000, 20_000, 50_000, 100_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_ECAL_Resolution_Gamma_9mm_Fe"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_Gamma_9mm_Fe"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_9mm_Fe[:,0])
ECAL.perform_fit(energies, plotting_directory)