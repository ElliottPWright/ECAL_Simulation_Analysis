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
            integrated_energy = np.sum(edeps)
            print(f"Integrated energy = {integrated_energy:.2f}")

            # All layers energy histogram
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(ak.to_numpy(ak.flatten(edeps)), bins = 100)
                plt.vlines(x=105.658, ymin=0, ymax=2*10**5, colors='red', linestyles='--', label='$\mu^{-}$ mass')
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.yscale("log")
                plt.legend()
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
                plt.show()
            

            # First layer energy histogram 
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(first_layer_edeps, bins = 100)
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.yscale("log")
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
                plt.show()

            # Last layer energy histogram 
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist(last_layer_edeps, bins = 100)
                plt.xlabel("Energy [MeV]")
                plt.ylabel("Counts")
                plt.yscale("log")
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_final_layer_energy.png"))
                plt.show()

            # All layers xy-hits histogram 
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(ak.to_numpy(ak.flatten(x_globals)), ak.to_numpy(ak.flatten(y_globals)), bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_final_layer_energy.png"))
                plt.show()

            # First xy-hits histogram 
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(first_layer_x_hits, first_layer_y_hits, bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
                plt.show() 

            # First layer xy-hits histogram 
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.hist2d(last_layer_x_hits, last_layer_y_hits, bins = 100, cmap ='viridis', cmin =1)
                plt.xlabel("x [mm]")
                plt.ylabel("y [mm]")
                plt.colorbar(label = 'Counts')
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
                plt.show()
            
            # Acceptance energy-based
            with plt.style.context(['science', 'no-latex']):
                plt.rcParams['figure.dpi'] = 200
                plt.plot(2160, acceptance_energy)
                plt.xlabel("Crack size [mm]")
                plt.ylabel("Acceptance")
                #plt.savefig(os.path.join(plotting_directory, "SplitCalGap_first_layer_energy.png"))
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


    
