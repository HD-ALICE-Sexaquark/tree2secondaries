#!/bin/bash

# `tree2secondaries/scripts/deploy.sh` -- Deploy files to either PI or GSI server.

set -euo pipefail

print_usage() { echo "usage: ./deploy.sh <pi|gsi>"; }

# check environment (1)
if [[ -z ${T2DS_LOCAL_PATH:-} ]]; then echo "error: missing env. var. T2DS_LOCAL_PATH"; exit 1; fi

# command-line arguments
if [[ $# -ne 1 ]]; then print_usage; exit 1; fi
server=$1
[[ " pi gsi " == *" ${server:-} "* ]] || { echo "error: invalid server"; exit 1; }

# check environment (2) based on server opt
server_upper=${server^^} # uppercase
host_var="${server_upper}_USER_AND_HOST"
path_var="T2DS_AT_${server_upper}_PATH"
[[ -z ${!host_var:-} ]] && { echo "error: missing env. var. ${host_var}"; exit 1; }
[[ -z ${!path_var:-} ]] && { echo "error: missing env. var. ${path_var}"; exit 1; }
t2ds_remote_user_and_host=${!host_var}
t2ds_remote_path=${!path_var}

# ensure remote directories exists
ssh "${t2ds_remote_user_and_host}" mkdir -p "${t2ds_remote_path}"

rsync -avzR --delete \
    "${T2DS_LOCAL_PATH}"/./common \
    "${T2DS_LOCAL_PATH}"/./configs \
    "${T2DS_LOCAL_PATH}"/./include \
    "${T2DS_LOCAL_PATH}"/./scripts \
    "${T2DS_LOCAL_PATH}"/./src \
    "${T2DS_LOCAL_PATH}"/./CMakeLists.txt \
    "${t2ds_remote_user_and_host}":"${t2ds_remote_path}"/

# create marker file
ssh "${t2ds_remote_user_and_host}" touch "${t2ds_remote_path}/CONTENT-DEPLOYED"
