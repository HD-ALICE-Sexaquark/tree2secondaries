#!/bin/bash

# !! Don't execute this script directly, it is meant to be used by `analysis_wrapper.sh` !!

#SBATCH --partition=main
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000

# check environment
if [[ -z ${T2S_OUTPUT_DIR} ]]; then echo "missing T2S_OUTPUT_DIR"; exit 1; fi
if [[ -z ${T2S_BIN} ]]; then echo "missing T2S_BIN"; exit 1; fi
# -- batch options
if [[ -z ${SIM_SET} ]]; then echo "missing SIM_SET"; exit 1; fi
if [[ -z ${MODE} ]]; then echo "missing MODE"; exit 1; fi
if [[ -z ${DATA_TYPE} ]]; then echo "missing DATA_TYPE"; exit 1; fi
if [[ -z ${PRODUCTION_PATH} ]]; then echo "missing PRODUCTION_PATH"; exit 1; fi
# -- per run number options
if [[ -z ${UNROLLED_RUNS} ]]; then echo "missing UNROLLED_RUNS"; exit 1; fi

# convert space-separated strings into arrays (NOTE: because Slurm)
all_runs=(${UNROLLED_RUNS})

# prepare executable parameters
run_number=${all_runs[${SLURM_ARRAY_TASK_ID}]}
production_name=$(basename "${PRODUCTION_PATH}")
input_file=${PRODUCTION_PATH}/${SIM_SET}/AnalysisResults_${run_number}.root
output_dir=${T2S_OUTPUT_DIR}/${production_name}/${SIM_SET}/${MODE^}
mkdir -pv "${output_dir}"
output_file=${output_dir}/${run_number}.root
reaction_channel=${SIM_SET::1}
sexa_mass=${SIM_SET:1}

t2s_command="${T2S_BIN} -i ${input_file} -o ${output_file} ${MODE} ${DATA_TYPE} -c ${reaction_channel} -m ${sexa_mass}"
echo "analysis_single_run @ farm-pi :: executing \"${t2s_command}\""
${t2s_command}
