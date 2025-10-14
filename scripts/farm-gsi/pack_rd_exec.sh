#!/bin/bash

# !! Don't execute this script directly, it is meant to be used by `pack_rd_wrapper.sh` !!

#SBATCH --partition=long
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000

# check environment
if [[ -z ${LUSTRE_HOME} ]]; then echo "missing LUSTRE_HOME"; exit 1; fi
# -- esd2vector
if [[ -z ${E2V_ROOT_DIR} ]]; then echo "missing E2V_ROOT_DIR"; exit 1; fi
if [[ -z ${E2V_OUTPUT_DIR} ]]; then echo "missing E2V_OUTPUT_DIR"; exit 1; fi
# -- tree2secondaries
if [[ -z ${T2S_OUTPUT_DIR} ]]; then echo "missing T2S_OUTPUT_DIR"; exit 1; fi
if [[ -z ${T2S_SLURM_DIR} ]]; then echo "missing T2S_SLURM_DIR"; exit 1; fi
if [[ -z ${T2S_BIN} ]]; then echo "missing T2S_BIN"; exit 1; fi
# -- main options
if [[ -z ${PRODUCTION_NAME} ]]; then echo "missing PRODUCTION_NAME"; exit 1; fi
if [[ -z ${MODE} ]]; then echo "missing MODE"; exit 1; fi
if [[ -z ${DATA_TYPE} ]]; then echo "missing DATA_TYPE"; exit 1; fi
# -- per run number options
if [[ -z ${UNROLLED_RUNS} ]]; then echo "missing UNROLLED_RUNS"; exit 1; fi
if [[ -z ${UNROLLED_DN} ]]; then echo "missing UNROLLED_DN"; exit 1; fi

# load software within container (hardcoded)
source /opt/root/v6-36-04/bin/thisroot.sh
export CMAKE_PREFIX_PATH=/opt/kfparticle/dev

# convert space-separated strings into arrays (NOTE: because Slurm)
all_runs=(${UNROLLED_RUNS})
all_dirs=(${UNROLLED_DN})

# extract info from job index
run_number=${all_runs[${SLURM_ARRAY_TASK_ID}]}
dir_number=${all_dirs[${SLURM_ARRAY_TASK_ID}]}

# prepare executable parameters
production_path=${E2V_OUTPUT_DIR}/${PRODUCTION_NAME}
input_files="${production_path}/${run_number}/AnalysisResults_${dir_number}.*.root"
output_dir=${T2S_OUTPUT_DIR}/${MODE}ed/${PRODUCTION_NAME}/${run_number}
mkdir -p "${output_dir}"
output_file="${output_dir}/${MODE^}ed_${dir_number}.root"

# log hack (https://unix.stackexchange.com/a/585453)
tmp_logfile=${T2S_SLURM_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log # = ${T2S_SLURM_DIR}/tmp/%A_%a.log
slurm_subdir=${T2S_SLURM_DIR}/${MODE}ed_${PRODUCTION_NAME}
mkdir -p "${slurm_subdir}"
ln -f "${tmp_logfile}" "${slurm_subdir}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

t2s_command="${T2S_BIN} -i ${input_files} -o ${output_file} ${MODE} ${DATA_TYPE}"
echo "executing \"${t2s_command}\""
${t2s_command}

rm "${tmp_logfile}" # end of log hack
