import os
import subprocess
import ROOT
import csv
import shutil
import numpy as np

# Function to update calo.cfg
def update_config(config_file, gap):
    with open(config_file, 'r') as file:
        lines = file.readlines()

    with open(config_file, 'w') as file:
        i = 0
        for line in lines:
            if 'airgap_mm' in line:
                file.write(f'airgap_mm = {gap}')
            else:
                file.write(line)


# -----------------------------
# Shell script executed on node
# -----------------------------




# Function to run the simulation
def run_simulation(executable):
    subprocess.run(executable)

def main():

    config_file = '/eos/user/e/elwright/GeoModelSplitCal/calo.cfg'
    executable = '/eos/user/e/elwright/GeoModelSplitCal/build/run_g4'

    gap = 1000
    #gaps = np.linspace(1000, 5000, 5)

    for energy in energies:
        print(f'Running simulation with gap: {gap} mm')
        update_config(config_file, gap)
        run_simulation(executable)


if __name__ == '__main__':
    main()
