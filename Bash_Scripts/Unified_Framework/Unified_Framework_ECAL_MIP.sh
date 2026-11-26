#!/bin/bash

export PATH=/cvmfs/ship-nightlies.cern.ch/pixi/0.72.0-rattlerfs-layers-d5a7568/bin:$PATH

cd ~/SHiP_Work/aegir

pixi cvmfs

cd ~/ECAL_Simulation_Analysis/Python

pixi run python Unified_Framework_Analysis.py
