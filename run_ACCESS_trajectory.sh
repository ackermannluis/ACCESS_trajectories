#!/bin/bash

#PBS -P k10
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l ncpus=37
#PBS -l mem=50gb
#PBS -l wd
#PBS -l storage=scratch/k10+gdata/k10+gdata/lb4
#PBS -l jobfs=2gb
#PBS -W umask=0022

set -eu

PYTHON_SCRIPT_NAME="ACCESS_trajectories_MP_V01.py"
trajectories_start_data_filename="ACCESS_trajectories_start_data.txt"
SCRIPT_ROOT=/g/data/k10/la6753/scripts


export PATH="${SCRIPT_ROOT}:${PATH}"

echo "starting to load build.env"
source ${SCRIPT_ROOT}/build.env
echo "build.env ran"

INPUT_DIR=/g/data/k10/la6753/scripts/$PYTHON_SCRIPT_NAME
cp ${INPUT_DIR} ${PBS_JOBFS}/
echo "copied script to node"

INPUT_DIR=/g/data/k10/la6753/scripts/$trajectories_start_data_filename
cp ${INPUT_DIR} ${PBS_JOBFS}/
echo "copied trajectories_start_data_filename to node"

INPUT_DIR=/g/data/k10/la6753/scripts/U_Analysis_main.py
cp ${INPUT_DIR} ${PBS_JOBFS}/

cd ${PBS_JOBFS}/


python3 $PYTHON_SCRIPT_NAME
echo "python script ran"

echo "finished."

