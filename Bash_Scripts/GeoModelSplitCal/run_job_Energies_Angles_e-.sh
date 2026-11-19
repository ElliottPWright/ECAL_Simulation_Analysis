#!/bin/bash

PROCESS=$1

# Energy array
ENERGIES=(
1000 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
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


cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_Energies_Angles_e_2/sim_mu-_${ENERGY}MeV_1000mm.root
fi

