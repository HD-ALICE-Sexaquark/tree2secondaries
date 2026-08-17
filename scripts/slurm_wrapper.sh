#!/bin/bash

# `tree2secondaries/scripts/slurm_wrapper.sh`
# ===========================================
# Slurm submission script for pack/search/verify, for both MC and real data.
#
# Expected input layout (both servers):
# - for data:
#   - for pack|verify:
#     PRODUCTION_PATH / <run number> / AnalysisResults_<file identifier>.root (multiple files are merged per chunk)
#   - for search:
#     PRODUCTION_PATH / PackedRNT_<run number>.root
# - for mc:
#   - for pack|verify:
#     PRODUCTION_PATH[_CHANNEL+MASS if any] / AnalysisResults_<run number>.root
#   - for search:
#     PRODUCTION_PATH[_CHANNEL+MASS if any] / PackedRNT_<run number>.root
#
# Output:
#   T2DS_ROOT_DIR / output / [packed,found,verified] / PRODUCTION_NAME[_CHANNEL+MASS if any] \
#                          / [Packed,Found,Verified]RNT_<run number><_chunk number if any>.root
#
# The submission writes a manifest, one line per array task, of the form:
#   <full t2ds command>
# `slurm_exec.sh` picks the line matching its array index and runs it verbatim.

set -euo pipefail
shopt -s nullglob

# server detection
export SERVER_NAME="pi"
if [[ ${HOSTNAME} != "alice-serv15" ]]; then
    SERVER_NAME="gsi"
fi

# hardcoded options
max_parallel_jobs=72 # pi only
chunk_size=25
slurm_time_merged="08:00:00" # pack|verify data: several files merged per task
slurm_time_single="01:00:00" # everything else: one file per task

print_usage() {
    echo "usage: ./slurm_wrapper.sh [-y] pack   <mc|data> <PRODUCTION PATH> <CHANNELS> <MASSES>"
    echo "       ./slurm_wrapper.sh [-y] search <mc|data> <PRODUCTION PATH> <CHANNELS> <MASSES>"
    echo "       ./slurm_wrapper.sh [-y] verify <mc|data> <PRODUCTION PATH> <CHANNELS> <MASSES>"
    echo "where: -y:              skip the confirmation of the executable's last modification time"
    echo "       PRODUCTION PATH: (see input layout description in script file (head -25))"
    echo "       CHANNELS:        comma-separated reaction channels (e.g. \"A,D\"), or \"\" for none"
    echo "                        - required different from \"\" for \`pack mc\`, \`search mc\`, \`search data\`"
    echo "                        - must be \"\" for \`pack data\` and \`verify data\`"
    echo "                        - can be different from \"\" for \`verify mc\` (for dirs selection)"
    echo "       MASSES:          comma-separated injected masses (e.g. \"1.73,1.8\"), or \"\" for none"
    echo "                        - not used at all in \`data\`"
    echo "                        - must be different from \"\" for \`pack mc\`"
    echo "                        - can be different from \"\" for \`verify mc\` (for dirs selection)"
}

# command-line arguments (1)
skip_bin_confirmation=0
while getopts ":y" opt; do
    case ${opt} in
        y) skip_bin_confirmation=1 ;;
        *) print_usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# command-line arguments (2)
if [[ $# -ne 5 ]]; then print_usage; exit 1; fi
mode=$1
if [[ ${mode} != "pack" && ${mode} != "search" && ${mode} != "verify" ]]; then print_usage; exit 1; fi
data_type=$2
if [[ ${data_type} != "mc" && ${data_type} != "data" ]]; then print_usage; exit 1; fi
# -- determine production name (NOTE: absolute, the manifest must not depend on the submission dir)
production_path=$(readlink -f "$3")
production_name=$(basename "${production_path}")
production_name=${production_name##*_} # remove prefix
# -- get reaction channels
IFS=',' read -ra reaction_channels <<< "$4"
if [[ ${#reaction_channels[@]} -eq 0 ]]; then
    reaction_channels=("NONE")
fi
# -- get injected masses
IFS=',' read -ra injected_masses <<< "$5"
if [[ ${#injected_masses[@]} -eq 0 ]]; then
    injected_masses=("NONE")
fi

# cross-check data type/mode arguments
if [[ ${mode} == "pack" && ${data_type} == "mc" && (${reaction_channels[0]} == "NONE" || ${injected_masses[0]} == "NONE") ]]; then
    echo "error: pack mc requires CHANNELS and MASSES"; print_usage; exit 1
elif [[ ${mode} == "search" && ${reaction_channels[0]} == "NONE" ]]; then
    echo "error: search mode requires CHANNELS"; print_usage; exit 1
fi

# determine data type/mode-related options
case ${mode} in
    pack)   stage="packed" ;;
    search) stage="found" ;;
    verify) stage="verified" ;;
esac
slurm_time=${slurm_time_single}
if [[ ${data_type} == "data" && ${mode} != "search" ]]; then slurm_time=${slurm_time_merged}; fi

# check environment (common to both servers)
if [[ -z ${T2DS_ROOT_DIR:-} ]]; then echo "error: missing env. var. T2DS_ROOT_DIR"; exit 1; fi
mkdir -p "${T2DS_ROOT_DIR}/output"
mkdir -p "${T2DS_ROOT_DIR}/slurm"
if [[ -z ${T2DS_BIN:-} ]]; then echo "error: missing env. var. T2DS_BIN"; exit 1; fi
if [[ ! -e ${T2DS_BIN} ]]; then echo "error: missing file ${T2DS_BIN}"; exit 1; fi

# check environment (gsi-specific)
if [[ ${SERVER_NAME} == "gsi" ]]; then
    if [[ -z ${LUSTRE_HOME:-} ]]; then echo "error: missing env. var. LUSTRE_HOME"; exit 1; fi
    if [[ -z ${T2DS_SINGULARITY_IMG:-} ]]; then echo "error: missing env. var. T2DS_SINGULARITY_IMG"; exit 1; fi
fi

# -- prepare the manifest (one line per array task)
export MANIFEST
MANIFEST="${T2DS_ROOT_DIR}/slurm/manifests/${mode}_${data_type}_${production_name}_$(date '+%Y%m%d-%H%M%S').txt"
mkdir -p "$(dirname "${MANIFEST}")"
: > "${MANIFEST}" # empty file if it already exists

# main loop
for reaction_channel in "${reaction_channels[@]}"; do
    for injected_mass in "${injected_masses[@]}"; do

        # -- define simulation set for sexa mc, regarless of data type and mode
        simset_suffix=""
        if [[ ${reaction_channel} != "NONE" && ${injected_mass} != "NONE" ]]; then
            simset_suffix="_${reaction_channel}${injected_mass}"
        fi

        # -- if simset is defined, define output suffix accordingly; if not, fill with reaction channel, if any
        output_suffix=""
        if [[ -n ${simset_suffix} ]]; then
            output_suffix="${simset_suffix}"
        elif [[ ${reaction_channel} != "NONE" ]]; then
            output_suffix="_${reaction_channel}"
        fi

        # -- subcommand options: `pack mc` takes -c and -m, both `search` take -c, the rest take none
        t2ds_options=""
        if [[ ${mode} == "pack" && ${data_type} == "mc" ]]; then
            t2ds_options="-c ${reaction_channel} -m ${injected_mass}"
        elif [[ ${mode} == "search" ]]; then
            t2ds_options="-c ${reaction_channel}"
        fi

        # -- define paths
        input_dir="${production_path}${simset_suffix}"
        output_dir="${T2DS_ROOT_DIR}/output/${stage}/${production_name}${output_suffix}"
        mkdir -p "${output_dir}"

        echo "slurm_wrapper @ ${HOSTNAME} :: adding files from ${input_dir}"

        if [[ ${data_type} == "data" && ${mode} != "search" ]]; then
            # -- one group per run-number dir, merged in chunks of at most `chunk_size` files
            for run_dir in "${input_dir}"/*/; do
                # (pack data, verify data)
                input_files=("${run_dir}"AnalysisResults_*.root)
                if [[ ${#input_files[@]} -eq 0 ]]; then continue; fi
                run_number=$(basename "${run_dir}")

                offset=0
                chunk=0
                while [[ ${offset} -lt ${#input_files[@]} ]]; do
                    output_id="${run_number}"
                    if [[ ${#input_files[@]} -gt ${chunk_size} ]]; then
                        output_id="${run_number}_${chunk}"
                    fi
                    {
                        printf '%q -i ' "${T2DS_BIN}"
                        printf '%q ' "${input_files[@]:${offset}:${chunk_size}}"
                        printf -- '-o %q %s %s %s\n' \
                            "${output_dir}/${stage^}RNT_${output_id}.root" "${mode}" "${data_type}" "${t2ds_options}"
                    } >> "${MANIFEST}"
                    offset=$((offset + chunk_size))
                    chunk=$((chunk + 1))
                done
            done
        else
            # -- one input file per task: mc holds one file per run number, and so do the packed files
            if [[ ${mode} == "search" ]]; then
                # (search mc, search data)
                input_files=("${input_dir}"/PackedRNT_*.root)
            else
                # (pack mc, verify mc)
                input_files=("${input_dir}"/AnalysisResults_*.root)
            fi

            for input_file in "${input_files[@]}"; do
                output_id=$(basename "${input_file}" .root)
                output_id=${output_id#*_} # remove prefix
                printf '%q -i %q -o %q %s %s %s\n' \
                    "${T2DS_BIN}" "${input_file}" \
                    "${output_dir}/${stage^}RNT_${output_id}.root" "${mode}" "${data_type}" "${t2ds_options}" \
                    >> "${MANIFEST}"
            done
        fi
    done
done

n_total_jobs=$(wc -l < "${MANIFEST}")
if [[ ${n_total_jobs} -eq 0 ]]; then
    echo "slurm_wrapper @ ${HOSTNAME} :: no input files found, nothing to submit"
    exit 1
fi

# -- print manifest path
echo "slurm_wrapper @ ${HOSTNAME} :: manifest: ${MANIFEST}"

# -- confirm state of executable
last_mod_bin=$(date -d @"$(stat -c %Y "${T2DS_BIN}")" '+%Y-%m-%d %H:%M:%S')
echo "slurm_wrapper @ ${HOSTNAME} :: executable ${T2DS_BIN} was last edited on ${last_mod_bin}."
if [[ ${skip_bin_confirmation} -eq 0 ]]; then
    read -p "slurm_wrapper @ ${HOSTNAME} :: do you want to continue? (y/n) " -r bin_confirmation
    if [[ ${bin_confirmation} != "y" ]]; then exit 1; fi
fi

array_max=$((n_total_jobs - 1))

# tmp dir for log hack
tmp_slurm_dir="${T2DS_ROOT_DIR}/slurm/tmp"
mkdir -p "${tmp_slurm_dir}"

exec_script="${T2DS_ROOT_DIR}/scripts/slurm_exec.sh"
if [[ ${SERVER_NAME} == "pi" ]]; then
    sbatch \
        --output="${tmp_slurm_dir}/%A_%a.log" \
        --array="0-${array_max}%${max_parallel_jobs}" \
        --time="${slurm_time}" \
        -- "${exec_script}" && \
    echo "slurm_wrapper @ ${HOSTNAME} :: a total of ${n_total_jobs} jobs have been submitted"
else
    sbatch \
        --output="${tmp_slurm_dir}/%A_%a.log" \
        --array="0-${array_max}" \
        --time="${slurm_time}" \
        --singularity-container="${T2DS_SINGULARITY_IMG}" \
        -- "${exec_script}" && \
    echo "slurm_wrapper @ ${HOSTNAME} :: a total of ${n_total_jobs} jobs have been submitted"
fi
