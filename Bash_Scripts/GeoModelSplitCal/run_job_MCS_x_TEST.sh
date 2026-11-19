#!/bin/bash

PROCESS=$1

# x-angle array for 0.5 degrees and 1 degrees
ANGLES=(
	872 1745
)

ANGLE=${ANGLES[$PROCESS]}

if [ $PROCESS -ge ${#ANGLES[@]} ]; then
    echo "Invalid process number: $PROCESS"
    exit 1
fi

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/

cd $SCRATCH/GeoModelSplitCal

CONFIG_FILE="run.cfg"
sed -i "s/^.*direction.*$/direction = ${ANGLE} 0 100000/" ${CONFIG_FILE}


cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 

if compgen -G "*.root" > /dev/null; then
 mkdir -p "/eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_MCS_x_10E5/"
 mv *.root "/eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_MCS_x_10E5/sim_mu-_10000MeV_${PROCESS}.root"
fi

