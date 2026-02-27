#! /usr/bin/env bash
# ------------------------------------------------------------------------------
# File:         compile_experiments.sh
# Description:  This script compiles the experiments required for the jobs.
# ------------------------------------------------------------------------------
# Determine whether the script is executed locally or on Lisa.
if [ $# -eq 0 ]; then
    echo "~ Running on Snellius.."

    module load 2024 gompi/2024a

    WDIR="${HOME}/ComplexityContesters/paper_PF_hardness"

    # Create tmp directories for results.
    mkdir -p "${TMPDIR}/${job}/results/"

    # Copy existing results if there are any.
    if [[ -d $WDIR/jobs/${job}/results/ ]]; then
        cp -r "${WDIR}/jobs/${job}/results" "${TMPDIR}/${job}/"
    fi
else
    echo "~ Running locally.."
    WDIR="."
fi

# Echo place where the results will finally be stored.
echo "Setting WDIR=$WDIR"

# Setup paths for compilatoin.
CFLAGS="-O3 -Wall -std=c++17"
PROSPR_PATH="${WDIR}/code/prospr/prospr/core/src"
EXP_PATH="${WDIR}/code/experiments"

# Remove old logs.
rm -rf build.log

# Compile dfs_bnb_mpi.
EXP_CODE="dfs_bnb_mpi"
echo -e "\nCompiling $EXP_CODE.." 2>&1 | tee -a build.log
mpic++ $CFLAGS -o "${EXP_CODE}" "${EXP_PATH}/${EXP_CODE}.cpp" \
    "${PROSPR_PATH}/utils.cpp" \
    "${PROSPR_PATH}/depth_first.cpp" "${PROSPR_PATH}/depth_first_bnb.cpp" \
    "${PROSPR_PATH}/protein.cpp" "${PROSPR_PATH}/amino_acid.cpp" \
    2>&1 | tee -a build.log

# Compile solve_remaining.
EXP_CODE="solve_remaining"
echo -e "\n\nCompiling $EXP_CODE.." 2>&1 | tee -a build.log
mpic++ $CFLAGS -o "${EXP_CODE}" "${EXP_PATH}/${EXP_CODE}.cpp" \
    "${PROSPR_PATH}/utils.cpp" \
    "${PROSPR_PATH}/depth_first.cpp" "${PROSPR_PATH}/depth_first_bnb.cpp" \
    "${PROSPR_PATH}/protein.cpp" "${PROSPR_PATH}/amino_acid.cpp" \
    2>&1 | tee -a build.log
