#!/bin/bash

PROCESS=$1

# Energy array
ENERGIES=(
10000 50000
)

ENERGY=${ENERGIES[$PROCESS]}

if [ $PROCESS -ge ${#ENERGIES[@]} ]; then
    echo "Invalid process number: $PROCESS"
    exit 1
fi

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/

cd $SCRATCH/GeoModelSplitCal

CONFIG_FILE="run.cfg"
sed -i "s/^.*energy_MeV.*$/energy_MeV = ${ENERGY}/" ${CONFIG_FILE}
sed -i "s/^.*particle.*$/particle = mu-/" ${CONFIG_FILE}
sed -i "s/^.*n_events.*$/n_events = 10000/" ${CONFIG_FILE}

# # Set the beam to be middle of one hextant
# sed -i "s/^.*position_mm.*$/position_mm = -1080 2160 90000/" ${CONFIG_FILE}

# Set the beam to be middle of 1 cm strip and half-ish in a 6 cm strip
sed -i "s/^.*position_mm.*$/position_mm = -1105 2185 90000/" ${CONFIG_FILE}


cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/experiment/ship/user/Mainz_Sim_Data/ewrightl/sim_runs/run_MIP_Calibration_9mm_Fe/sim_mu-_${ENERGY}MeV_1000mm.root
fi

