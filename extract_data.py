#!/usr/bin/env python3
import glob
import re
import numpy as np
import uproot
import awkward as ak

EOS_PATH = "/eos/user/m/mahajanv/GeoModelSplitCal/sim_results/sim_e-_10k_*MeV.root"

def extract_energy(filename):
    match = re.search(r'(\d+)MeV', filename)
    return int(match.group(1)) if match else None

file_list = sorted(glob.glob(EOS_PATH))

for filepath in file_list:
    energy_mev = extract_energy(filepath)
    if energy_mev is None:
        continue

    print(f"Opening {filepath}...")
    with uproot.open(filepath) as f:
        tree = f["calo_events"]
        # Reads 'total_energy' branch into numpy array
        edep_hits = tree["edep"].array(library="ak")
        edep_arr = ak.to_numpy(ak.sum(edep_hits, axis=-1))
    
    out_file = f"data/edep_{energy_mev}MeV.npy"
    np.save(out_file, edep_arr)
    print(f"Saved {len(edep_arr)} entries -> {out_file}")
