#!/bin/bash

PROCESS=$1

# Gap array
GAPS=(
0
250 500
)

GAP=${GAPS[$PROCESS]}

if [ $PROCESS -ge ${#GAPS[@]} ]; then
    echo "Invalid process number: $PROCESS"
    exit 1
fi

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/

CONFIG_FILE="calo.cfg"
sed -i "s/^.*airgap_mm.*$/airgap_mm = ${GAP}/" ${CONFIG_FILE}


cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_TEST/sim_mu-_10000MeV_${GAP}mm.root
fi

