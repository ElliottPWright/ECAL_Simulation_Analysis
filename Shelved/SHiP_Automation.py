from __future__ import print_function

import argparse
import os
import subprocess


# -----------------------------
# Shell script executed on node
# -----------------------------
SHELL_SCRIPT = """
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


cd $SCRATCH/GeoModelSplitCal/build

chmod +x run_g4

./run_g4


if compgen -G "*.root" > /dev/null; then
mv *.root /eos/user/e/elwright/EPWL_ECAL_Simulations/sim_runs/run_MIP_Calibration/sim_mu-_${ENERGY}MeV_1000mm.root
fi

"""


# -----------------------------
# Submit file (NO EOS paths here except output storage)
# -----------------------------
SUBMIT_FILE = """
# HTCondor submit file for CERN LXPLUS
# Job: SHiP_Calo_test

universe   = vanilla
executable = run_job_MIP_Calibration.sh

arguments  = $(Process)

# Resource requirements
request_cpus   = 1
request_memory = 8GB
request_disk   = 4GB


log    = logs/calo_production.$(Cluster).$(Process).log
output = logs/calo_production.$(Cluster).$(Process).out
error  = logs/calo_production.$(Cluster).$(Process).err

# Job duration - automatically set based on runtime_hours=2.5h
+MaxRuntime             = 9000
+JobFlavour             = "workday"


# Notify user on completion (optional)
# notification            = Complete
# notify_user             = elwright@cern.ch

queue 2
"""


# -----------------------------
# Build submit + script
# -----------------------------


def create_submit_script(output_dir, njobs, tag=""):

    os.makedirs("logs", exist_ok=True)
    out_dir = os.path.join(output_dir, tag)
    script_name = f".run_simulation_{abs(hash(output_dir))}.sh"

    submit_text = SUBMIT_FILE.format(
        script_id=abs(hash(output_dir)),
        out_dir=os.path.abspath(out_dir),
        njobs=njobs,
    )

    return submit_text, script_name


# -----------------------------
# Main submit function
# -----------------------------
def submit():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-n", type=int, required=True)
    parser.add_argument("-t", "--tag", default="run")

    args = parser.parse_args()

    # MUST run from AFS, not EOS
    if os.getcwd().startswith("/eos"):
        raise RuntimeError("Do NOT submit from EOS. Run from AFS or home directory.")

    sub_file, script_name = create_submit_script(
        njobs=args.n,
        tag=args.tag,
        output_dir=args.output
    )

    with open("submit.sub", "w") as f:
        f.write(sub_file)

    with open(script_name, "w") as f:
        f.write(SHELL_SCRIPT)

    os.chmod(script_name, 0o755)

    print("Submitting jobs...")
    subprocess.call(["condor_submit", "submit.sub"])
    print("Done")


if __name__ == "__main__":
    submit()
                       
