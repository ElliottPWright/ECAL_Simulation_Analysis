import os
import uproot
import numpy as np
import awkward as ak
from pathlib import Path
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
import scienceplots

class SimFileReader:
    "This class opens and reads a simulation ROOT file."

    def __init__(self, file_path: str):
        self.path = Path(file_path)
        self.file = uproot.open(file_path)

    def list_trees(self):
       return [k for k in self.file.keys()]

    def get_branches(self, branches):
        "Get the relevant TTree branches from the ROOT file."

        trees = self.list_trees()
        tree = self.file[trees[0]]  # Assuming the first tree is the one we want
        df = tree.arrays(branches, library="ak")

        data = [df[branch] for branch in branches]

        return data 

    def Diagnostic_Plotting(self, file_path):
            # First of all, need to extract the relevant branches

            Sim = SimFileReader(os.path.join(file_path, "sim_mu-_10000MeV_2160mm.root"))
            
            branches = ["edep", "x_global", "y_global", "z_global", "layer"]
            edeps, x_globals, y_globals, z_globals, layers = Sim.get_branches(branches)

            # Now concatenate all events using awkward

            # edeps = ak.to_numpy(ak.flatten(edeps))
            # x_globals = ak.to_numpy(ak.flatten(x_globals))
            # y_globals = ak.to_numpy(ak.flatten(y_globals))
            # z_globals = ak.to_numpy(ak.flatten(z_globals))
            # layers = ak.to_numpy(ak.flatten(layers))

            first_layer_edeps = [edeps[i][layers[i] == 1] for i in range(len(edeps))]
            last_layer_edeps = [edeps[i][layers[i]==np.max(layers[i])] for i in range(len(edeps))]

            ## First layer energy histogram 
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(first_layer_edeps, bins = 100)
                plt.xlabel("Energy [GeV]")
                plt.ylabel("Counts")
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
                plt.show()

            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(last_layer_edeps, bins = 100)
                plt.xlabel("Energy [GeV]")
                plt.ylabel("Counts")
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_final_layer_energy.png"))
                plt.show()

    def close(self):
        self.file.close()


class ProcessHits(SimFileReader):
    "Basic hits processing."
    def __init__(self):
        self.first_layer_x_hits = self.x_globals[self.layers == 1]
        self.last_layer_x_hits = self.x_globals[self.layers == np.max(self.layers)]

        self.acceptance_energy = np.sum(self.last_layer_x_hits)/np.sum(self.first_layer_x_hits)
        self.acceptance_number = len(self.last_layer_x_hits)/len(self.first_layer_x_hits)


    
