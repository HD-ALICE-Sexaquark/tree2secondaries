#!/bin/bash

# `tree2secondaries/scripts/slurm_wrapper.sh`
# ===========================================
# Slurm submission script for pack/search/verify, for both MC and real data (DATA_TYPE=mc|data).
#
# Expected input layout (both servers):
#   - for data:
#     PRODUCTION_PATH / <run number> / AnalysisResults_<file identifier>.root (pack, verify)
#   - for mc:
#     PRODUCTION_PATH[_CHANNEL+MASS or nothing] / AnalysisResults_<file identifier>.root (pack, verify)
#     PRODUCTION_PATH[_CHANNEL+MASS or nothing] / PackedRNT_<file identifier>.root       (search)
#
# Output:
#   T2DS_ROOT_DIR / output / [packed,found,verified] / PRODUCTION_NAME[_CHANNEL+MASS or nothing] / [Packed,Found,Verified]RNT_<file identifier>.root

set -euo pipefail
shopt -s nullglob

# server detection
export SERVER_NAME="pi"
if [[ ${HOSTNAME} != "alice-serv15" ]]; then
    SERVER_NAME="gsi"
fi

# hardcoded options
max_parallel_jobs=72 # pi only

print_usage() {
    echo "usage: ./slurm_wrapper.sh pack   <mc|data> <PRODUCTION PATH> <CHANNELS> <MASSES>"
    echo "       ./slurm_wrapper.sh search <mc|data> <PRODUCTION PATH> <CHANNELS> <MASSES>"
    echo "       ./slurm_wrapper.sh verify <mc|data> <PRODUCTION PATH> <CHANNELS> <MASSES>"
    echo "where: CHANNELS: comma-separated reaction channels (e.g. \"A,D\"), or \"\" for none"
    echo "       MASSES:   comma-separated injected masses (e.g. \"1.73,1.8\"), or \"\" for none -- only used when DATA_TYPE=mc"
    echo "       (for DATA_TYPE=mc, CHANNELS and MASSES must either both be empty or both be non-empty)"
}

# command-line arguments
if [[ $# -ne 5 ]]; then print_usage; exit 1; fi
export MODE=$1
if [[ ${MODE} != "pack" && ${MODE} != "search" && ${MODE} != "verify" ]]; then print_usage; exit 1; fi
export DATA_TYPE=$2
if [[ ${DATA_TYPE} != "mc" && ${DATA_TYPE} != "data" ]]; then print_usage; exit 1; fi
export PRODUCTION_PATH=$3
export PRODUCTION_NAME
PRODUCTION_NAME=$(basename "${PRODUCTION_PATH}")
PRODUCTION_NAME=${PRODUCTION_NAME##*_} # remove long prefix

IFS=',' read -ra reaction_channels <<< "$4"
if [[ ${DATA_TYPE} == "mc" ]]; then
    IFS=',' read -ra injected_masses <<< "$5"
    if [[ ${#reaction_channels[@]} -eq 0 && ${#injected_masses[@]} -eq 0 ]]; then
        reaction_channels=("NONE")
        injected_masses=("NONE")
    elif [[ ${#reaction_channels[@]} -eq 0 || ${#injected_masses[@]} -eq 0 ]]; then
        print_usage; exit 1
    fi
else
    injected_masses=("NONE") # data has no injected masses
    if [[ ${#reaction_channels[@]} -eq 0 ]]; then reaction_channels=("NONE"); fi
fi

# check environment (common to both servers)
if [[ -z ${T2DS_ROOT_DIR:-} ]]; then echo "error: missing env. var. T2DS_ROOT_DIR"; exit 1; fi
mkdir -p "${T2DS_ROOT_DIR}/output"
mkdir -p "${T2DS_ROOT_DIR}/slurm"
if [[ -z ${T2DS_BIN:-} ]]; then echo "error: missing env. var. T2DS_BIN"; exit 1; fi
if [[ ! -e ${T2DS_BIN} ]]; then echo "error: missing file ${T2DS_BIN}"; exit 1; fi

# check environment (server-specific)
if [[ ${SERVER_NAME} == "gsi" ]]; then
    if [[ -z ${LUSTRE_HOME:-} ]]; then echo "error: missing env. var. LUSTRE_HOME"; exit 1; fi
    if [[ -z ${T2DS_SINGULARITY_IMG:-} ]]; then echo "error: missing env. var. T2DS_SINGULARITY_IMG"; exit 1; fi # PLACEHOLDER: path unknown yet
fi

# -- confirm state of executable
last_mod_bin=$(date -d @"$(stat -c %Y "${T2DS_BIN}")" '+%Y-%m-%d %H:%M:%S')
echo "slurm_wrapper @ ${HOSTNAME} :: executable ${T2DS_BIN} was last edited on ${last_mod_bin}."
read -p "slurm_wrapper @ ${HOSTNAME} :: do you want to continue? (y/n) " -r bin_confirmation
if [[ ${bin_confirmation} != "y" ]]; then exit 1; fi

# define strings (NOTE: not arrays, because Slurm)
export CHANNELS_STR
export MASSES_STR
export FILE_IDS_STR
n_total_jobs=0

# main loop
for reaction_channel in "${reaction_channels[@]}"; do
    for injected_mass in "${injected_masses[@]}"; do
        if [[ ${reaction_channel} == "NONE" ]]; then
            sim_set=""
            echo "adding files directly from ${PRODUCTION_PATH}"
        elif [[ ${DATA_TYPE} == "mc" ]]; then
            sim_set="${reaction_channel}${injected_mass}"
            echo "adding sim_set ${sim_set}"
        else
            sim_set="${reaction_channel}"
            echo "adding sim_set ${sim_set}"
        fi

        if [[ ${MODE} == "search" ]]; then
            input_glob="${PRODUCTION_PATH}${sim_set:+_${sim_set}}/PackedRNT_*.root"
        else
            # mode == "pack" or "verify"
            if [[ ${DATA_TYPE} == "mc" ]]; then
                input_glob="${PRODUCTION_PATH}${sim_set:+_${sim_set}}/AnalysisResults_*.root"
            else
                # data type == "data"
                input_glob="${PRODUCTION_PATH}/*/AnalysisResults_*.root"
            fi
        fi

        for fid_file in ${input_glob}; do

            # file id:
            # - for mc: <run number>
            # - for data, first: <output dir number>
            file_identifier=$(basename "${fid_file}" .root)
            file_identifier=${file_identifier##*_} # remove long prefix
            if [[ ${DATA_TYPE} == "data" ]]; then
                rd_run_number=$(basename "$(dirname "${fid_file}")")
                file_identifier=${rd_run_number}_${file_identifier} # data, second: <run number>_<output dir number>
            fi

            CHANNELS_STR+="${reaction_channel} "
            MASSES_STR+="${injected_mass} "
            FILE_IDS_STR+="${file_identifier} "

            n_total_jobs=$((n_total_jobs + 1))
        done
    done
done

array_max=$((n_total_jobs - 1))

# tmp dir for log hack
tmp_slurm_dir="${T2DS_ROOT_DIR}/slurm/tmp"
mkdir -p "${tmp_slurm_dir}"

exec_script="${T2DS_ROOT_DIR}/scripts/slurm_exec.sh"
if [[ ${n_total_jobs} -gt 0 ]]; then
    if [[ ${SERVER_NAME} == "pi" ]]; then
        sbatch \
            --output="${tmp_slurm_dir}/%A_%a.log" \
            --array="0-${array_max}%${max_parallel_jobs}" \
            -- "${exec_script}" && \
        echo "slurm_wrapper @ ${HOSTNAME} :: a total of ${n_total_jobs} jobs have been submitted"
    else
        sbatch \
            --output="${tmp_slurm_dir}/%A_%a.log" \
            --array="0-${array_max}" \
            --singularity-container="${T2DS_SINGULARITY_IMG}" \
            -- "${exec_script}" && \
        echo "slurm_wrapper @ ${HOSTNAME} :: a total of ${n_total_jobs} jobs have been submitted"
    fi
fi
