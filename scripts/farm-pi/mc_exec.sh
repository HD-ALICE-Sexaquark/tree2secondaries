#!/bin/bash

# `tree2secondaries/scripts/farm-pi/mc_exec.sh`
# =============================================
# NOTE: Don't execute this script directly, it is meant to be used by `tree2secondaries/scripts/farm-pi/mc_wrapper.sh`

#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

# strict mode
set -euo pipefail

# check environment
if [[ -z ${T2S_OUTPUT_DIR:-} ]]; then echo "error: missing env. var. T2S_OUTPUT_DIR"; exit 1; fi
if [[ -z ${T2S_BIN:-} ]]; then echo "error: missing env. var. T2S_BIN"; exit 1; fi
if [[ -z ${T2S_SLURM_DIR:-} ]]; then echo "error: missing env. var. T2S_SLURM_DIR"; exit 1; fi
# -- batch options
if [[ -z ${MODE:-} ]]; then echo "error: missing env. var. MODE"; exit 1; fi
if [[ -z ${DATA_TYPE:-} ]]; then echo "error: missing env. var. DATA_TYPE"; exit 1; fi
if [[ -z ${PRODUCTION_NAME:-} ]]; then echo "error: missing env. var. PRODUCTION_NAME"; exit 1; fi
if [[ ${MODE} != "search" && -z ${PRODUCTION_PATH:-} ]]; then echo "error: missing env. var. PRODUCTION_PATH"; exit 1; fi
# -- per run number options
if [[ -z ${CHANNELS_STR:-} ]]; then echo "error: missing env. var. CHANNELS_STR"; exit 1; fi
if [[ -z ${MASSES_STR:-} ]]; then echo "error: missing env. var. MASSES_STR"; exit 1; fi
if [[ -z ${RUN_NUMBERS_STR:-} ]]; then echo "error: missing env. var. RUN_NUMBERS_STR"; exit 1; fi

# convert space-separated strings into arrays (NOTE: because Slurm)
read -ra CHANNELS_ARR <<< "${CHANNELS_STR}"
read -ra MASSES_ARR <<< "${MASSES_STR}"
read -ra RUN_NUMBERS_ARR <<< "${RUN_NUMBERS_STR}"

# extract info from job index
reaction_channel=${CHANNELS_ARR[${SLURM_ARRAY_TASK_ID}]}
sexa_mass=${MASSES_ARR[${SLURM_ARRAY_TASK_ID}]}
run_number=${RUN_NUMBERS_ARR[${SLURM_ARRAY_TASK_ID}]}

# prepare executable parameters
sim_set="${reaction_channel}${sexa_mass}"
case ${MODE} in
    pack)
        stage="packed"
        input_file=${PRODUCTION_PATH}/${sim_set}/AnalysisResults_${run_number}.root
        t2s_options="-c ${reaction_channel} -m ${sexa_mass}"
        ;;
    search)
        stage="found"
        input_file=${T2S_OUTPUT_DIR}/packed/${PRODUCTION_NAME}_${sim_set}/Packed_${run_number}.root
        t2s_options="-c ${reaction_channel}"
        ;;
    verify)
        stage="verified"
        input_file=${PRODUCTION_PATH}/${sim_set}/AnalysisResults_${run_number}.root
        t2s_options=""
        ;;
    *)
        echo "error: unknown MODE ${MODE}"; exit 1
        ;;
esac
output_dir=${T2S_OUTPUT_DIR}/${stage}/${PRODUCTION_NAME}_${sim_set}
mkdir -p "${output_dir}"
output_file=${output_dir}/${stage^}_${run_number}.root

# log hack (https://unix.stackexchange.com/a/585453)
tmp_logfile=${T2S_SLURM_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log # = ${T2S_SLURM_DIR}/tmp/%A_%a.log
slurm_subdir=${T2S_SLURM_DIR}/${stage}_${PRODUCTION_NAME}_${sim_set}
mkdir -p "${slurm_subdir}"
ln -f "${tmp_logfile}" "${slurm_subdir}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

t2s_command="${T2S_BIN} -i ${input_file} -o ${output_file} ${MODE} ${DATA_TYPE} ${t2s_options}"
echo "executing \"${t2s_command}\""
${t2s_command}

rm "${tmp_logfile}" # end of log hack
