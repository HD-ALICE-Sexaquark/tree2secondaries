#!/bin/bash

# `tree2secondaries/scripts/farm-pi/mc_wrapper.sh`
# ================================================
# Script to submit batch Slurm jobs for the three operations: pack, search and verify.
# - `pack` and `verify` work if the MC output of `esd2vector` has the following dir structure:
#   [production path] / [simulation set] / AnalysisResults_[run number].root
# - `search` works if the packed MC of `tree2secondaries` has the following dir structure:
#   T2S_OUTPUT_DIR / packed / [production name]_[simulation set] / Packed_[run number].root
# The output will be:
#   T2S_OUTPUT_DIR / [packed,found,verified] / [production name]_[simulation set] / [Packed,Found,Verified]_[run number].root

# strict mode
set -euo pipefail
shopt -s nullglob

# hardcoded options
max_parallel_jobs=32
reaction_channels=("A") # ("A" "D")
injected_masses=("1.8") # ("1.73" "1.8" "1.87" "1.94" "2.01")
export DATA_TYPE="mc"

# function
print_usage() {
    echo "usage: ./mc_wrapper.sh pack <PRODUCTION PATH>"
    echo "       ./mc_wrapper.sh search <PRODUCTION NAME> (e.g. 23l1a3, 23l1b3)"
    echo "       ./mc_wrapper.sh verify <PRODUCTION PATH>"
}

# command-line arguments
if [[ $# -ne 2 ]]; then print_usage; exit 1; fi
export MODE=$1
if [[ ${MODE} != "pack" && ${MODE} != "search" && ${MODE} != "verify" ]]; then print_usage; exit 1; fi
export PRODUCTION_PATH=""
export PRODUCTION_NAME=""
if [[ ${MODE} == "search" ]]; then
    PRODUCTION_NAME=$2
else
    PRODUCTION_PATH=$2
    PRODUCTION_NAME=$(basename "${PRODUCTION_PATH}")
    PRODUCTION_NAME=${PRODUCTION_NAME##*_} # remove long prefix
fi

# check environment
if [[ -z ${T2S_OUTPUT_DIR:-} ]]; then echo "error: missing env. var. T2S_OUTPUT_DIR"; exit 1; fi
mkdir -p "${T2S_OUTPUT_DIR}"
if [[ -z ${T2S_SLURM_DIR:-} ]]; then echo "error: missing env. var. T2S_SLURM_DIR"; exit 1; fi
mkdir -p "${T2S_SLURM_DIR}"
# -- confirm state of executable
if [[ -z ${T2S_BIN:-} ]]; then echo "error: missing env. var. T2S_BIN"; exit 1; fi
if [[ ! -e ${T2S_BIN} ]]; then echo "error: missing file ${T2S_BIN}"; exit 1; fi
last_mod_bin=$(date -d @"$(stat -c %Y "${T2S_BIN}")" '+%Y-%m-%d %H:%M:%S')
echo "mc_wrapper @ ${HOSTNAME} :: executable ${T2S_BIN} was last edited on ${last_mod_bin}."
read -p "mc_wrapper @ ${HOSTNAME} :: do you want to continue? (y/n) " -r bin_confirmation
if [[ ${bin_confirmation} != "y" ]]; then exit 1; fi

# define strings (NOTE: not arrays, because Slurm)
export CHANNELS_STR
export MASSES_STR
export RUN_NUMBERS_STR
n_total_jobs=0

# main loop
for reaction_channel in "${reaction_channels[@]}"; do
    for injected_mass in "${injected_masses[@]}"; do
        sim_set="${reaction_channel}${injected_mass}"
        echo "adding sim_set ${sim_set}"

        if [[ ${MODE} == "search" ]]; then
            input_glob="${T2S_OUTPUT_DIR}/packed/${PRODUCTION_NAME}_${sim_set}/Packed_*.root"
        else
            input_glob="${PRODUCTION_PATH}/${sim_set}/AnalysisResults_*.root"
        fi

        for rn_file in ${input_glob}; do

            run_number=$(basename "${rn_file}" .root)
            run_number=${run_number##*_} # remove long prefix
            echo "  adding rn ${run_number}"

            CHANNELS_STR+="${reaction_channel} "
            MASSES_STR+="${injected_mass} "
            RUN_NUMBERS_STR+="${run_number} "

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
        --output="${tmp_slurm_dir}/%A_%a.log" \
        --array="0-${array_max}%${max_parallel_jobs}" \
        -- mc_exec.sh && \
    echo "mc_wrapper @ ${HOSTNAME} :: a total of ${n_total_jobs} jobs have been submitted"
fi
