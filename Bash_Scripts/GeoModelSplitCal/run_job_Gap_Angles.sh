#!/bin/bash

PROCESS=$1

# Gap array
GAPS=(
0
250 500 750 1000 1250 1500 1750 2000
2250 2500 2750 3000 3250 3500 3750 4000
4250 4500 4750 5000 5250 5500 5750 6000
6250 6500 6750 7000 7250 7500 7750 8000
8250 8500 8750 9000 9250 9500 9750 10000
)
# Layer codes: 1=WidePVT, 3=ThinPS, 5 HPL,  7=Lead, 8 split, maybe there is a better way to input these codes? (Below are the original organised into groups)
# layers = 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2, 5,6, 8, 
#          7,3,7,4, 7,1,7,2, 5,6, 7,3,7,4 ,7,1,7,2,7,3,7,4 ,7,1,7,2, 5,6, 7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4

# Pb is used as the absorber in the ECAL
# lead_thickness_mm = 3

GAP=${GAPS[$PROCESS]}

if [ $PROCESS -ge ${#GAPS[@]} ]; then
    echo "Invalid process number: $PROCESS"
    exit 1
fi

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/

cd $SCRATCH/GeoModelSplitCal

CONFIG_FILE="calo.cfg"
sed -i "s/^.*airgap_mm.*$/airgap_mm = ${GAP}/" ${CONFIG_FILE}


cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4 


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_Gap_Angles/sim_mu-_10000MeV_${GAP}mm.root
fi

