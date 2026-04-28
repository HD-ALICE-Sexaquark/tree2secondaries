#!/bin/bash

# `tree2secondaries/scripts/download_files.sh`
# ============================================
# Download output files from the preceding program `esd2vector`.

# strict mode
set -euo pipefail

# functions
print_usage() {
    echo "USAGE:"
    echo "  ./download_files.sh --server <pi|gsi> --cycle <cycle> --n <n>";
    echo "where:"
    echo "  <cycle> -- the name of the dirs between \"output/\" and the AnalysisResults_<rn>.root files"
    echo "  <n>     -- max number of files to download";
    echo "EXAMPLES:"
    echo "  ./download_files.sh --server pi --cycle local_mc_23l1a3/A1.8 --n 5"
}
process_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --server)  server="$2"; shift 2 ;;
            --cycle)   cycle="$2"; shift 2 ;;
            --n)       max_n_files="$2"; shift 2 ;;
            *)         echo "error: unrecognized argument: $1"; print_usage; exit 1 ;;
        esac
    done
}

# check environment (1)
if [[ -z ${T2S_LOCAL_PATH:-} ]]; then echo "error: missing env. var. T2S_LOCAL_PATH"; exit 1; fi

# command-line arguments
if [[ $# -lt 6 ]]; then print_usage; exit 1; fi
process_args "$@"

# validate input args
[[ " pi gsi " == *" ${server:-} "* ]] || { echo "error: invalid server"; exit 1; }
[[ ${max_n_files} =~ ^[0-9]+$ ]] || { echo "error: --n must be a positive integer"; exit 1; }

# check environment (2) based on server opt
server_upper=${server^^}
host_var="${server_upper}_USER_AND_HOST"
path_var="E2V_AT_${server_upper}_PATH"
[[ -z ${!host_var:-} ]] && { echo "error: missing env. var. ${host_var}"; exit 1; }
[[ -z ${!path_var:-} ]] && { echo "error: missing env. var. ${path_var}"; exit 1; }
remote_user_and_host=${!host_var}
remote_top_dir=${!path_var}

# main #

local_files_dir=${T2S_LOCAL_PATH}/files/
mkdir -p "${local_files_dir}"

mapfile -t list_remote_files < <(ssh -q "${remote_user_and_host}" "ls '${remote_top_dir}/output/${cycle}'")

files_counter=0
for file in "${list_remote_files[@]}"; do
    if [[ ${files_counter} -ge ${max_n_files} ]]; then break; fi
    rsync -avzR \
        "${remote_user_and_host}":"${remote_top_dir}/output/./${cycle}/${file}" \
        "${local_files_dir}"
    files_counter=$((files_counter + 1))
done
