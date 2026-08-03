#!/bin/bash

PROCESS=$1

# Gap array
CRACKS=(
2160 2170 2180 2210
)

ENERGIES=(
10000 20000 30000 50000
)



CRACK=${CRACKS[$(( PROCESS % 4 ))]}
ENERGY=${ENERGIES[$(( PROCESS / 4 ))]}


source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/

cd $SCRATCH/GeoModelSplitCal

CONFIG_FILE="calo.cfg"
sed -i "s/^.*module_pitch_y_mm .*$/module_pitch_y_mm  = ${CRACK}/" ${CONFIG_FILE}

RUN_FILE="run.cfg"
sed -i "s/^.*energy_MeV .*$/energy_MeV = ${ENERGY}/" ${RUN_FILE}

cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_Cracks_Uniform/sim_mu-_${ENERGY}MeV_${CRACK}mm.root
fi

