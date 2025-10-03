#!/bin/bash

# Script to submit batch Slurm jobs.
# It works if the MC output of `esd2vector` has the following dir structure:
# [production dir]
# ├── [simulation set]
# │   ├── AnalysisResults_[run number 1].root
# │   └── AnalysisResults_[run number 2].root

function print_usage() { echo "Usage: pack_mc_wrapper.sh [production dir]"; }

# hardcoded options
max_parallel_jobs=32
reaction_channels=("A" "D")
injected_masses=("1.73" "1.8" "1.87" "1.94" "2.01")
export MODE="pack"
export DATA_TYPE="mc"

# check environment
if [[ -z ${T2S_OUTPUT_DIR} ]]; then echo "missing T2S_OUTPUT_DIR"; exit 1; fi
mkdir -p "${T2S_OUTPUT_DIR}"
if [[ -z ${T2S_SLURM_DIR} ]]; then echo "missing T2S_SLURM_DIR"; exit 1; fi
mkdir -p "${T2S_SLURM_DIR}"
# -- confirm state of executable
if [[ -z ${T2S_BIN} ]]; then echo "missing T2S_BIN"; exit 1; fi
if [[ ! -e ${T2S_BIN} ]]; then echo "missing file ${T2S_BIN}"; exit 1; fi
echo "pack_mc_wrapper @ farm-pi :: Executable ${T2S_BIN} was last edited on $(stat -c %y "${T2S_BIN}"). Do you want to continue? (y/n)"
read -n 1 -p "y" answer
if [[ ${answer} != "y" ]]; then exit 1; fi

# command-line arguments
if [[ $# -ne 1 ]]; then print_usage; exit 1; fi
export PRODUCTION_PATH=$1

# define strings (NOTE: not arrays, because Slurm)
export UNROLLED_CHANNELS=""
export UNROLLED_MASSES=""
export UNROLLED_RUNS=""
n_total_jobs=0

# main loop
for reaction_channel in "${reaction_channels[@]}"; do
    for injected_mass in "${injected_masses[@]}"; do
        sim_set="${reaction_channel}${injected_mass}"
        for rn_file in "${PRODUCTION_PATH}/${sim_set}"/AnalysisResults_*.root; do

            run_number=${rn_file/.root/} # remove .root
            run_number=${run_number##*_} # remove long prefix

            UNROLLED_CHANNELS+="${reaction_channel} "
            UNROLLED_MASSES+="${injected_mass} "
            UNROLLED_RUNS+="${run_number} "

            n_total_jobs=$((n_total_jobs + 1))
        done
    done
done

array_max=$((n_total_jobs - 1))

# tmp dir for log hack
tmp_slurm_dir="${T2S_SLURM_DIR}/tmp"
mkdir -p "${tmp_slurm_dir}"

if [[ ${n_total_jobs} -gt 0 ]]; then
    sbatch \
        --output="${tmp_slurm_dir}"/%A_%a.log \
        --array="0-${array_max}%${max_parallel_jobs}" \
        -- pack_mc_exec.sh && \
    echo "pack_mc_wrapper @ farm-pi :: a total of ${n_total_jobs} jobs have been submitted"
fi
