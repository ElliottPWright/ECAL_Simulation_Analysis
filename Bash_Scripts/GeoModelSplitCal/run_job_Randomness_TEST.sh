#!/bin/bash

PROCESS=$1

X=$(( (RANDOM % 2000) - 1000 - 1600 ))
Y=$(( (RANDOM % 3000) - 1500 - 1600 ))
Z=33000

DX=$(( (RANDOM % 4000) - $X - 2000 ))
DY=$(( (RANDOM % 6000) - $Y - 3000 ))
DZ=50000

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/

cd $SCRATCH/GeoModelSplitCal

RUN_FILE="run.cfg"

sed -i "s/^.*position_mm.*$/position_mm = ${X} ${Y} ${Z}/" ${RUN_FILE}
sed -i "s/^.*direction.*$/direction = ${DX} ${DY} ${DZ}/" ${RUN_FILE}

cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_Randomness_TEST/sim_mu-_10000MeV_${PROCESS}mm.root
fi

