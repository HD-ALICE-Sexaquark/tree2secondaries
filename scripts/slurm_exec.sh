#!/bin/bash

# `tree2secondaries/scripts/slurm_exec.sh`
# =========================================
# NOTE: Don't execute this script directly, it is meant to be used by `tree2secondaries/scripts/slurm_wrapper.sh`

#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

set -euo pipefail

# check environment (common to both servers)
if [[ -z ${T2DS_ROOT_DIR:-} ]]; then echo "error: missing env. var. T2DS_ROOT_DIR"; exit 1; fi
if [[ -z ${T2DS_BIN:-} ]]; then echo "error: missing env. var. T2DS_BIN"; exit 1; fi
if [[ -z ${SERVER_NAME:-} ]]; then echo "error: missing env. var. SERVER_NAME"; exit 1; fi
# -- batch options
if [[ -z ${MODE:-} ]]; then echo "error: missing env. var. MODE"; exit 1; fi
if [[ -z ${DATA_TYPE:-} ]]; then echo "error: missing env. var. DATA_TYPE"; exit 1; fi
if [[ -z ${PRODUCTION_NAME:-} ]]; then echo "error: missing env. var. PRODUCTION_NAME"; exit 1; fi
if [[ -z ${PRODUCTION_PATH:-} ]]; then echo "error: missing env. var. PRODUCTION_PATH"; exit 1; fi
# -- per-job options
if [[ -z ${CHANNELS_STR:-} ]]; then echo "error: missing env. var. CHANNELS_STR"; exit 1; fi
if [[ -z ${MASSES_STR:-} ]]; then echo "error: missing env. var. MASSES_STR"; exit 1; fi
if [[ -z ${FILE_IDS_STR:-} ]]; then echo "error: missing env. var. FILE_IDS_STR"; exit 1; fi

# environment (server-specific)
if [[ ${SERVER_NAME} == "gsi" ]]; then
    if [[ -z ${LUSTRE_HOME:-} ]]; then echo "error: missing env. var. LUSTRE_HOME"; exit 1; fi
    # load software within container (NOTE: hardcoded to the current image)
    source /opt/root/v6-36-04/bin/thisroot.sh
    export CMAKE_PREFIX_PATH=/opt/kfparticle/dev
fi

# convert space-separated strings into arrays (NOTE: because Slurm)
read -ra CHANNELS_ARR <<< "${CHANNELS_STR}"
read -ra MASSES_ARR <<< "${MASSES_STR}"
read -ra FILE_IDS_ARR <<< "${FILE_IDS_STR}"

# extract info from job index
reaction_channel=${CHANNELS_ARR[${SLURM_ARRAY_TASK_ID}]}
sexa_mass=${MASSES_ARR[${SLURM_ARRAY_TASK_ID}]}
file_identifier=${FILE_IDS_ARR[${SLURM_ARRAY_TASK_ID}]}

# prepare executable parameters
if [[ ${reaction_channel} == "NONE" ]]; then
    sim_set=""
elif [[ ${DATA_TYPE} == "mc" ]]; then
    sim_set="${reaction_channel}${sexa_mass}"
else
    sim_set="${reaction_channel}"
fi
name_suffix=${sim_set:+_${sim_set}}
case ${MODE} in
    pack)
        stage="packed"
        input_file=${PRODUCTION_PATH}${name_suffix}/AnalysisResults_${file_identifier}.root
        t2ds_options=""
        if [[ ${reaction_channel} != "NONE" ]]; then
            t2ds_options="-c ${reaction_channel}"
            if [[ ${DATA_TYPE} == "mc" ]]; then t2ds_options+=" -m ${sexa_mass}"; fi
        fi
        ;;
    search)
        stage="found"
        input_file=${PRODUCTION_PATH}${name_suffix}/PackedRNT_${file_identifier}.root
        t2ds_options=""
        if [[ ${reaction_channel} != "NONE" ]]; then t2ds_options="-c ${reaction_channel}"; fi
        ;;
    verify)
        stage="verified"
        if [[ ${DATA_TYPE} == "mc" ]]; then
            input_file=${PRODUCTION_PATH}${name_suffix}/AnalysisResults_${file_identifier}.root
        else
            # data type == "data"
            run_number=${file_identifier%_*}
            out_dir_number=${file_identifier#*_}
            input_file=${PRODUCTION_PATH}/${run_number}/AnalysisResults_${out_dir_number}.root
        fi
        t2ds_options=""
        ;;
    *)
        echo "slurm_exec @ ${HOSTNAME} :: ERROR :: unknown MODE ${MODE}"; exit 1
        ;;
esac
output_dir=${T2DS_ROOT_DIR}/output/${stage}/${PRODUCTION_NAME}${name_suffix}
mkdir -p "${output_dir}"
output_file=${output_dir}/${stage^}RNT_${file_identifier}.root

# log hack (https://unix.stackexchange.com/a/585453)
tmp_logfile=${T2DS_ROOT_DIR}/slurm/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log # = ${T2DS_ROOT_DIR}/slurm/tmp/%A_%a.log
slurm_subdir=${T2DS_ROOT_DIR}/slurm/${stage}_${PRODUCTION_NAME}${name_suffix}
mkdir -p "${slurm_subdir}"
ln -f "${tmp_logfile}" "${slurm_subdir}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

t2ds_command="${T2DS_BIN} -i ${input_file} -o ${output_file} ${MODE} ${DATA_TYPE} ${t2ds_options}"
echo "executing \"${t2ds_command}\""
${t2ds_command}

rm "${tmp_logfile}" # end of log hack
