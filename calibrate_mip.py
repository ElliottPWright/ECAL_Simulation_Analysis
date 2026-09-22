#!/usr/bin/env python3

import os
import numpy as np
import uproot
from collections import defaultdict

# ============================================================
# Input
# ============================================================

FILE = "/eos/user/m/mahajanv/GeoModelSplitCal/sim_results/sim_mu-_10k_5000MeV.root"

# Output directory
os.makedirs("data", exist_ok=True)


# ============================================================
# Read ROOT file
# ============================================================

print(f"Opening {FILE}...")

with uproot.open(FILE) as f:
    tree = f["calo_events"]

    edep = tree["edep"].array(library="np")
    typ = tree["type"].array(library="np")
    sec = tree["section"].array(library="np")
    lay = tree["layer"].array(library="np")
    vol = tree["vol"].array(library="np")


print(f"Number of events: {len(edep)}")


# ============================================================
# Collect strip energies
#
# type 1,2 = 6 cm WidePVT strips
# type 3,4 = 1 cm ThinPS strips
#
# A strip is uniquely identified by:
# (type, section, layer, volume)
#
# Hits belonging to the same strip are summed.
# ============================================================

strip_1cm = []
strip_6cm = []

for event_edep, event_type, event_section, event_layer, event_vol in zip(
    edep, typ, sec, lay, vol
):

    # Sum energy deposits belonging to the same strip
    strips = defaultdict(float)

    for energy, strip_type, section, layer, volume in zip(
        event_edep,
        event_type,
        event_section,
        event_layer,
        event_vol
    ):

        key = (
            int(strip_type),
            int(section),
            int(layer),
            int(volume)
        )

        strips[key] += float(energy)

    # Store the total energy deposited in each strip
    for key, energy in strips.items():

        strip_type = key[0]

        # 1 cm strips
        if strip_type in (3, 4):
            strip_1cm.append(energy)

        # 6 cm strips
        elif strip_type in (1, 2):
            strip_6cm.append(energy)


# ============================================================
# Convert to NumPy arrays
# ============================================================

strip_1cm = np.asarray(strip_1cm, dtype=float)
strip_6cm = np.asarray(strip_6cm, dtype=float)


# ============================================================
# Print basic information
# ============================================================

print()
print("======================================")
print("Muon MIP calibration data")
print("======================================")

print(f"1 cm strip deposits : {len(strip_1cm)}")
print(f"6 cm strip deposits : {len(strip_6cm)}")

if len(strip_1cm) > 0:
    print()
    print("1 cm:")
    print(f"  min    = {np.min(strip_1cm):.6f} MeV")
    print(f"  max    = {np.max(strip_1cm):.6f} MeV")
    print(f"  mean   = {np.mean(strip_1cm):.6f} MeV")
    print(f"  median = {np.median(strip_1cm):.6f} MeV")

if len(strip_6cm) > 0:
    print()
    print("6 cm:")
    print(f"  min    = {np.min(strip_6cm):.6f} MeV")
    print(f"  max    = {np.max(strip_6cm):.6f} MeV")
    print(f"  mean   = {np.mean(strip_6cm):.6f} MeV")
    print(f"  median = {np.median(strip_6cm):.6f} MeV")


# ============================================================
# Save data
# ============================================================

file_1cm = "data/muon_1cm_strip_energy.npy"
file_6cm = "data/muon_6cm_strip_energy.npy"

np.save(file_1cm, strip_1cm)
np.save(file_6cm, strip_6cm)

print()
print("Saved:")
print(f"  {file_1cm}")
print(f"  {file_6cm}")
print()
print("Done.")
