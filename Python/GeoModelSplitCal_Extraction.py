import os
import uproot
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from GeoModelSplitCal_Extraction_Classes import SimFileReader, ProcessHits

def main_exec():
    file_path =  r"C:\Users\ewrightl\cernbox\EPWL_ECAL_Simulations\sim_runs\run_Cracks_Gaussian"
    Sim = SimFileReader(os.path.join(file_path, "sim_mu-_10000MeV_2160mm.root"))

    branches = ["edep", "x_global", "y_global", "z_global", "type"]
    edeps, x_globals, y_globals, z_globals, types = Sim.get_branches(branches)
    # print(f"Branches extracted: {branches}")
    # print(f"Number of events: {len(edeps)}")
    # print(f"Number of hits in first event: {len(edeps[1])}")
    # print(f"First event energy deposits: {edeps[1]}")
    # print(f"First event x hits: {x_globals[1]}")
    # print(f"First event y hits: {y_globals[1]}")
    # print(f"First event z hits: {z_globals[1]}")
    # print(f"First event layer types: {types[1]}")

    print(x_globals[:][types == 1])

    Sim.Diagnostic_Plotting(file_path)

if __name__ == "__main__":
    main_exec()
