#!/bin/bash

# !! Don't execute this script directly, it is meant to be used by `search_mc_wrapper.sh` !!

#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

# check environment
if [[ -z ${T2S_OUTPUT_DIR} ]]; then echo "missing T2S_OUTPUT_DIR"; exit 1; fi
if [[ -z ${T2S_BIN} ]]; then echo "missing T2S_BIN"; exit 1; fi
if [[ -z ${T2S_SLURM_DIR} ]]; then echo "missing T2S_SLURM_DIR"; exit 1; fi
# -- batch options
if [[ -z ${MODE} ]]; then echo "missing MODE"; exit 1; fi
if [[ -z ${DATA_TYPE} ]]; then echo "missing DATA_TYPE"; exit 1; fi
if [[ -z ${PRODUCTION_NAME} ]]; then echo "missing PRODUCTION_NAME"; exit 1; fi
# -- per run number options
if [[ -z ${UNROLLED_CHANNELS} ]]; then echo "missing UNROLLED_CHANNELS"; exit 1; fi
if [[ -z ${UNROLLED_MASSES} ]]; then echo "missing UNROLLED_MASSES"; exit 1; fi
if [[ -z ${UNROLLED_RUNS} ]]; then echo "missing UNROLLED_RUNS"; exit 1; fi

# convert space-separated strings into arrays (NOTE: because Slurm)
all_channels=(${UNROLLED_CHANNELS})
all_masses=(${UNROLLED_MASSES})
all_runs=(${UNROLLED_RUNS})

# extract info from job index
reaction_channel=${all_channels[${SLURM_ARRAY_TASK_ID}]}
sexa_mass=${all_masses[${SLURM_ARRAY_TASK_ID}]}
run_number=${all_runs[${SLURM_ARRAY_TASK_ID}]}

# prepare executable parameters
sim_set="${reaction_channel}${sexa_mass}"
input_dir=${T2S_OUTPUT_DIR}/packed/${PRODUCTION_NAME}_${sim_set}
input_file=${input_dir}/Packed_${run_number}.root
mode_str=; if [[ ${MODE} == "pack" ]]; then mode_str="packed"; else mode_str="found"; fi
output_dir=${T2S_OUTPUT_DIR}/${mode_str}/${PRODUCTION_NAME}_${sim_set}
mkdir -p "${output_dir}"
output_file=${output_dir}/${mode_str^}_${run_number}.root

# log hack (https://unix.stackexchange.com/a/585453)
tmp_logfile=${T2S_SLURM_DIR}/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log # = ${T2S_SLURM_DIR}/tmp/%A_%a.log
slurm_subdir=${T2S_SLURM_DIR}/${mode_str}_${PRODUCTION_NAME}_${sim_set}
mkdir -p "${slurm_subdir}"
ln -f "${tmp_logfile}" "${slurm_subdir}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

t2s_command="${T2S_BIN} -i ${input_file} -o ${output_file} ${MODE} ${DATA_TYPE} -c ${reaction_channel} -m ${sexa_mass}"
echo "executing \"${t2s_command}\""
${t2s_command}

rm "${tmp_logfile}" # end of log hack
