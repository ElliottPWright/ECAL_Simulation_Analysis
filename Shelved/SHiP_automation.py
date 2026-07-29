import os
import subprocess
import ROOT
import csv
import shutil
import numpy as np

# Function to update calo.cfg
def update_config(config_file, energy):
    with open(config_file, 'r') as file:
        lines = file.readlines()

    with open(config_file, 'w') as file:
        i = 0
        for line in lines:
            if 'energy_MeV' in line:
                file.write(f'energy_MeV = {energy}')
            else:
                file.write(line)

# Function to run the simulation
def run_simulation(executable):
    subprocess.run(executable)

def main():

    config_file = '/eos/user/e/elwright/GeoModelSplitCal/calo.cfg'
    executable = '/eos/user/e/elwright/GeoModelSplitCal/build/run_g4'


    energies = np.linspace(1000, 5000, 5)

    for energy in energies:
        print(f'Running simulation with muon energy: {energy} MeV')
        update_config(config_file, energy)
        run_simulation(executable)


if __name__ == '__main__':
    main()
