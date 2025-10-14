#!/bin/bash

# Script to submit batch Slurm jobs.
# It works if the Real Data (RD) output of `esd2vector` has the following dir structure:
#
# E2V_OUTPUT_DIR
# └── [production name]
#     └── [run number]
#         ├── AnalysisResults_[dir number A1].[dir number B1].root
#         └── AnalysisResults_[dir number A1].[dir number B2].root
#
# The output will be:
#
# T2S_OUTPUT_DIR
# └── packed
#     └── [production name]
#         └── [run number]
#             ├── Packed_[dir number A1].root
#             └── Packed_[dir number A2].root

function print_usage() { echo "Usage: pack_rd_wrapper.sh [LHC15o,LHC18q,LHC18r]"; }

# command-line arguments
if [[ $# -ne 1 ]]; then print_usage; exit 1; fi
export PRODUCTION_NAME=$1

# hardcoded options
max_parallel_jobs=500
export MODE="pack"
export DATA_TYPE="data"

# check environment
if [[ -z ${LUSTRE_HOME} ]]; then echo "missing LUSTRE_HOME"; exit 1; fi
# -- esd2vector
if [[ -z ${E2V_ROOT_DIR} ]]; then echo "missing E2V_ROOT_DIR"; exit 1; fi
if [[ -z ${E2V_OUTPUT_DIR} ]]; then echo "missing E2V_OUTPUT_DIR"; exit 1; fi
# -- tree2secondaries
if [[ -z ${T2S_OUTPUT_DIR} ]]; then echo "missing T2S_OUTPUT_DIR"; exit 1; fi
mkdir -p "${T2S_OUTPUT_DIR}"
if [[ -z ${T2S_SLURM_DIR} ]]; then echo "missing T2S_SLURM_DIR"; exit 1; fi
mkdir -p "${T2S_SLURM_DIR}"
# -- confirm state of executable
if [[ -z ${T2S_BIN} ]]; then echo "missing T2S_BIN"; exit 1; fi
if [[ ! -e ${T2S_BIN} ]]; then echo "missing file ${T2S_BIN}"; exit 1; fi
last_mod_bin=$(date -d @"$(stat -c %Y "${T2S_BIN}")" '+%Y-%m-%d %H:%M:%S')
echo "pack_mc_wrapper @ farm-pi :: Executable ${T2S_BIN} was last edited on ${last_mod_bin}."
read -p "pack_mc_wrapper @ farm-pi :: Do you want to continue? (y/n) " -r bin_confirmation
if [[ ${bin_confirmation} != "y" ]]; then exit 1; fi

# determine year and pass number
export YEAR_2DIG=${PRODUCTION_NAME:3:2}
export PASS_NUMBER=2
if [[ ${YEAR_2DIG} -eq 18 ]]; then PASS_NUMBER=3; fi

# define strings (NOTE: not arrays, because Slurm)
export UNROLLED_RUNS=""
export UNROLLED_DN=""
n_total_jobs=0

# loop over run numbers (NOTE: from text file, because Lustre)
run_numbers_file=${E2V_ROOT_DIR}/doc/${PRODUCTION_NAME}_pass${PASS_NUMBER}_rn.txt
while read -r line1; do

    run_number=${line1}
    # n_jobs_per_rn=0 # DEBUG

    # loop over dir numbers (NOTE: from text file, because Lustre)
    dir_numbers_file=${E2V_ROOT_DIR}/doc/dir_numbers/${PRODUCTION_NAME}/${run_number}.txt
    current_dn=0
    while read -r line2; do
        dn=${line2/%.*} # remove long suffix
        if [[ "${dn}" == "${current_dn}" ]]; then continue; fi
        current_dn=${dn}
        # echo "pack_rd_wrapper @ farm-gsi :: adding jobs for RN ${run_number} and all DN ${current_dn}.*" # DEBUG

        UNROLLED_RUNS+="${run_number} "
        UNROLLED_DN+="${current_dn} "

        n_total_jobs=$((n_total_jobs + 1))
        # n_jobs_per_rn=$((n_jobs_per_rn + 1)) # DEBUG
    done < "${dir_numbers_file}"

    # echo "pack_rd_wrapper @ farm-gsi :: RN ${run_number} = ${n_jobs_per_rn} remaining jobs" # DEBUG
done < "${run_numbers_file}"

array_max=$((n_total_jobs - 1))
mkdir -p "${T2S_SLURM_DIR}/tmp"

if [[ ${n_total_jobs} -gt 0 ]]; then
    sbatch \
        --singularity-container="${LUSTRE_HOME}/containers/root+kf.sif" \
        --output="${T2S_SLURM_DIR}/tmp/%A_%a.log" \
        --array="0-${array_max}%${max_parallel_jobs}" \
        -- pack_rd_exec.sh && \
    echo "pack_rd_wrapper @ farm-gsi :: a total of ${n_total_jobs} jobs have been submitted"
fi
