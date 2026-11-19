#!/bin/bash

PROCESS=$1

# LAYERS[0] = Original design
# LAYERS[1] = Remove 1 ECAL layer before the gap (-7,1,7,2,7,3,7,4,) and increase Pb thickness by 52mm/36 = 1.44 mm
# LAYERS[2] = Remove 2 ECAL layers (1 before, 1 after gap, 2x -7,1,7,2,7,3,7,4,) and increase Pb thickness by 104mm/32 = 3.25 mm
# LAYERS[3] = Remove 3 ECAL layers (1 before, 2 after gap, 3x -7,1,7,2,7,3,7,4,) and increase Pb thickness by 156mm/28 = 5.57 mm
# LAYERS[4] = Remove 4 ECAL layers (2 before, 2 after gap, 4x -7,1,7,2,7,3,7,4,) and increase Pb thickness by 208mm/24 = 8.67 mm

# Gap array
LAYERS=(
    "7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,5,6,8,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4"
    "7,1,7,2,7,3,7,4,7,1,7,2,5,6,8,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4"
    "7,1,7,2,7,3,7,4,7,1,7,2,5,6,8,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4"
    "7,1,7,2,7,3,7,4,7,1,7,2,5,6,8,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4"
    "7,1,7,2,5,6,8,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,5,6,7,3,7,4,7,1,7,2,7,3,7,4,7,1,7,2,7,3,7,4"

)

ENERGIES=(1000 2000 4000 6000 10000)

THICKNESSES=(3 4.44 7.25 8.57 11.67)



EIDX=$(( PROCESS / 5 ))
IDX=$(( PROCESS % 5 ))


THICKNESS=${THICKNESSES[$IDX]}

LAYER=${LAYERS[$IDX]}

ENERGY=${ENERGIES[$EIDX]}

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh

SCRATCH=$_CONDOR_SCRATCH_DIR

mkdir -p $SCRATCH

cp -r /eos/user/e/elwright/GeoModelSplitCal $SCRATCH/


cd $SCRATCH/GeoModelSplitCal

CONFIG_FILE="calo.cfg"
RUN_FILE="run.cfg"

for i in {1}
do
    sed -i "s/^.*layers.*$/layers = ${LAYER}/" ${CONFIG_FILE}
done

sed -i "s/^.*lead_thickness_mm.*$/lead_thickness_mm = ${THICKNESS}/" ${CONFIG_FILE}

sed -i "s/^.*energy_MeV.*$/energy_MeV = ${ENERGY}/" ${RUN_FILE}

cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_Layers_Energy_2/sim_mu-_${ENERGY}MeV_${IDX}.root
fi


# Layer codes: 1=WidePVT, 3=ThinPS, 5 HPL,  7=Lead, 8 split, maybe there is a better way to input these codes? (Below are the original organised into groups)
# layers = 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2, 5,6, 8, 
#          7,3,7,4, 7,1,7,2, 5,6, 7,3,7,4 ,7,1,7,2,7,3,7,4 ,7,1,7,2, 5,6, 7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4, 7,1,7,2,7,3,7,4
