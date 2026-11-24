from Unified_Framework_Extraction_Classes import SimFileReader, ProcessHits, MIP_Calibration, ECAL_Resolution
import ROOT
import os
import types
import sys
import numpy as np

calibration = [2.30, 0.03]
input_directory = r"/afs/cern.ch/user/e/elwright/SHiP_Work/outputs/tutorial_muon_beam.root"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/Unified_Framework"


# MIP calibration

Sim = SimFileReader(input_directory)
MIP = MIP_Calibration(input_directory, plotting_directory)

calibration_3mm_Pb = MIP.perform_fit()

print()
print(".........................................")
print("Printing Results of MIP calibration")
print(calibration_3mm_Pb)


# ECAL e- Resolution

energies = np.array([10_000])

file_path =  r"/eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_ECAL_Resolution"
plotting_directory = r"/afs/cern.ch/user/e/elwright/ECAL_Simulation_Analysis/Figures/ECAL_Resolution_3mm_Pb"

ECAL = ECAL_Resolution(file_path, plotting_directory, calibration_3mm_Pb[:,0])

results = ECAL.perform_fit(energies, plotting_directory)



## Below is just for testing

# r = ROOT.RNTupleReader.Open("events", os.path.join(directory, "muon_gun.root")) # or "straw_tubes_hits", ...
# r.PrintInfo(); r.Show(0)

_alive = []

## Here follows procedural code that will be implemented into
## the class structure in due course 

# Reading one file

path = os.path.join(directory, "tutorial_muon_beam.root")

f = SimFileReader(path)
f.view("event_number", "std::uint64_t")
f.view("sim_hits", "std::vector<SHiP::SimHit>")
f.view("sim_particles", "std::vector<SHiP::SimParticle>")

z_positions = []
edeps = []

for i in f.entries():
    hits = f["sim_hits"](i)

    for h in hits:
        z_position = h.position[2]
        z_positions.append(z_position)

        if z_position < 500:
            edep = h.energy_deposit
            edeps.append(edep)
    

# ## Looking at the energy deposited in the ECAL

# counts, edges = np.histogram(edeps, bins=100)
# h = ROOT.TH1D("h_edep", "", len(counts), edges)
# for j, c in enumerate(counts, start=1):
#     h.SetBinContent(j, c)

# canvas = ROOT.TCanvas("c_edep_fit", "edep_fit", 900, 700)
# canvas.SetLeftMargin(0.15)

# h.GetXaxis().SetTitle("Energy deposited per strip [GeV]")
# h.GetYaxis().SetTitle("Entries")

# # Change hist axis

# h.GetXaxis().SetRangeUser(0, 30)
# h.Draw()

# canvas.SaveAs("MIP_Calibration.png")

# ## Looking at the z-pos of the hits

# z_positions = np.asarray(z_positions).flatten()


# counts, edges = np.histogram(z_positions, bins=100)
# h = ROOT.TH1D("h_edep", "", len(counts), edges)
# for j, c in enumerate(counts, start=1):
#     h.SetBinContent(j, c)

# canvas = ROOT.TCanvas("c_z_positions", "z_positions", 900, 700)
# canvas.SetLeftMargin(0.15)

# h.GetXaxis().SetTitle("z hit position [cm]")
# h.GetYaxis().SetTitle("Entries")

# # Change hist axis

# h.Draw()

# canvas.SaveAs("z_positions.png")