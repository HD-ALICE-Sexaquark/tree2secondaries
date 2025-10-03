#!/bin/bash

function print_usage() {
    echo "Usage: analysis_wrapper.sh [-s] [-m] [-d] [production path]"
    echo ""
    echo "Options:"
    echo "      -s <reaction_channel><sexa_mass> -m <mode> -d <data_type>"
    echo ""
    echo "Where:"
    echo "      reaction_channel : A, D"
    echo "      sexa_mass        : 1.73, 1.8, 1.87, 1.94, 2.01"
    echo "      mode             : pack, search"
    echo "      data_type        : mc, data"
    echo ""
    echo "Example:"
    echo "      analysis_wrapper.sh -s A1.8 -m pack -d mc [production path]"
}

# check environment
if [[ -z ${T2S_ROOT_DIR} ]]; then echo "missing T2S_ROOT_DIR"; exit 1; fi
if [[ -z ${T2S_BIN} ]]; then echo "missing T2S_BIN"; exit 1; fi
if [[ -z ${T2S_OUTPUT_DIR} ]]; then echo "missing T2S_OUTPUT_DIR"; exit 1; fi
mkdir -pv "${T2S_OUTPUT_DIR}"
if [[ -z ${T2S_SLURM_DIR} ]]; then echo "missing T2S_SLURM_DIR"; exit 1; fi
mkdir -pv "${T2S_SLURM_DIR}"

# command-line arguments
# -- flags
export SIM_SET=
export MODE=
export DATA_TYPE=
while getopts s:m:d: flag; do
    case "${flag}" in
        s) SIM_SET=${OPTARG};;
        m) MODE=${OPTARG};;
        d) DATA_TYPE=${OPTARG};;
        *) print_usage; exit 1;;
    esac
done
# -- positional options
shift $((OPTIND - 1))
export PRODUCTION_PATH=$1

# hardcoded options
export MAX_PARALLEL_JOBS=24

# define strings (NOTE: not arrays, because Slurm)
export UNROLLED_RUNS=""
n_total_jobs=0

# loop over run numbers
for rn_file in "${PRODUCTION_PATH}"/"${SIM_SET}"/AnalysisResults_*.root; do

    run_number=${rn_file/.root/} # remove .root
    run_number=${run_number##*_} # remove long prefix

    UNROLLED_RUNS+="${run_number} "

    n_total_jobs=$((n_total_jobs + 1))
done

array_max=$((n_total_jobs - 1))

sub_slurm_dir="${T2S_SLURM_DIR}/${PRODUCTION_PATH}/${SIM_SET}"
mkdir -pv "${sub_slurm_dir}"

if [[ ${n_total_jobs} -gt 0 ]]; then
    sbatch \
        --output="${sub_slurm_dir}"/%A_%a.log \
        --array="0-${array_max}%${MAX_PARALLEL_JOBS}" \
        -- analysis_single_run.sh && \
    echo "analysis_wrapper @ farm-pi :: a total of ${n_total_jobs} jobs have been submitted"
fi
