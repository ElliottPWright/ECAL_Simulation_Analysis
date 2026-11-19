# Code based on Matei's Event class in the ShipSoft framework
# booklet. Adjusted for the Mainz ECAL use case.

import ROOT
import os
import types
import sys
import numpy as np
#import matplotlib.pyplot as plt


directory = r"/afs/cern.ch/user/e/elwright/SHiP_Work/outputs"

# r = ROOT.RNTupleReader.Open("events", os.path.join(directory, "muon_gun.root")) # or "straw_tubes_hits", ...
# r.PrintInfo(); r.Show(0)

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
    

## Looking at the energy deposited in the ECAL

counts, edges = np.histogram(edeps, bins=100)
h = ROOT.TH1D("h_edep", "", len(counts), edges)
for j, c in enumerate(counts, start=1):
    h.SetBinContent(j, c)

canvas = ROOT.TCanvas("c_edep_fit", "edep_fit", 900, 700)
canvas.SetLeftMargin(0.15)

h.GetXaxis().SetTitle("Energy deposited per strip [GeV]")
h.GetYaxis().SetTitle("Entries")

# Change hist axis

h.GetXaxis().SetRangeUser(0, 30)
h.Draw()

canvas.SaveAs("MIP_Calibration.png")

## Looking at the z-pos of the hits

z_positions = np.asarray(z_positions).flatten()


counts, edges = np.histogram(z_positions, bins=100)
h = ROOT.TH1D("h_edep", "", len(counts), edges)
for j, c in enumerate(counts, start=1):
    h.SetBinContent(j, c)

canvas = ROOT.TCanvas("c_z_positions", "z_positions", 900, 700)
canvas.SetLeftMargin(0.15)

h.GetXaxis().SetTitle("z hit position [cm]")
h.GetYaxis().SetTitle("Entries")

# Change hist axis

h.Draw()

canvas.SaveAs("z_positions.png")