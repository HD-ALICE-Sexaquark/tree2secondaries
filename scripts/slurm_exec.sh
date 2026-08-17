#!/bin/bash

# `tree2secondaries/scripts/slurm_exec.sh`
# =========================================
# NOTE: Don't execute this script directly, it is meant to be used by `tree2secondaries/scripts/slurm_wrapper.sh`
#
# One array task = one line of MANIFEST:
#   <full t2ds command, already built by the wrapper>

#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

set -euo pipefail

# check environment
if [[ -z ${T2DS_ROOT_DIR:-} ]]; then echo "error: missing env. var. T2DS_ROOT_DIR"; exit 1; fi
if [[ -z ${MANIFEST:-} ]]; then echo "error: missing env. var. MANIFEST"; exit 1; fi
if [[ -z ${T2DS_BIN:-} ]]; then echo "error: missing env. var. T2DS_BIN"; exit 1; fi

# prepare environment
# -- on gsi, the job runs within the container, but `sbatch` still exports the login node's environment:
#    drop the ROOT-related vars, as `LD_LIBRARY_PATH` outranks the executable's `RUNPATH`
if [[ ${SERVER_NAME:-} == "gsi" ]]; then unset LD_LIBRARY_PATH ROOTSYS PYTHONPATH; fi
# -- make `libPOD.so` and its rootmap reachable regardless of the working directory
t2ds_bin_dir=$(dirname "${T2DS_BIN}")
export LD_LIBRARY_PATH="${t2ds_bin_dir}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

# extract this task's command from the manifest (NOTE: only this task's line, never the whole file)
manifest_line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${MANIFEST}")
if [[ -z ${manifest_line} ]]; then
    echo "slurm_exec @ ${HOSTNAME} :: ERROR :: no line $((SLURM_ARRAY_TASK_ID + 1)) in ${MANIFEST}"; exit 1
fi

# log hack (https://unix.stackexchange.com/a/585453)
tmp_logfile=${T2DS_ROOT_DIR}/slurm/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log # = ${T2DS_ROOT_DIR}/slurm/tmp/%A_%a.log
slurm_subdir="${T2DS_ROOT_DIR}/slurm/$(basename "${MANIFEST}" .txt)"
mkdir -p "${slurm_subdir}"
ln -f "${tmp_logfile}" "${slurm_subdir}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

echo "slurm_exec @ ${HOSTNAME} :: executing line $((SLURM_ARRAY_TASK_ID + 1)) of ${MANIFEST}"
eval "${manifest_line}"

rm "${tmp_logfile}" # end of log hack; if something above crashed, the non-deletion of this file is a feature
