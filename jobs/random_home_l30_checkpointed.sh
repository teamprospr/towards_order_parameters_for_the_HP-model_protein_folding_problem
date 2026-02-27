#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# This script folds the proteins of length 30 for the random dataset.
# Pass an argument, doesn't matter what, to trigger local execution.
# The default is Snellius execution.
# With 128 tasks and 1 node, this script takes around ?? minutes to run.
# ------------------------------------------------------------------------------
#SBATCH --job-name=random_l30
#SBATCH --output=slurm/%x-%j.out
#SBATCH --error=slurm/%x-%j.err
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --ear=on
#SBATCH --ear-policy=monitoring
#SBATCH --ear-verbose=1
#SBATCH --mail-user=okke.van.eck@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
# Allow job to be requeued (re-uses output through checkpointing).
# Send SIGINT 10 seconds before end of job.
#SBATCH --signal=INT@10                    
#SBATCH --requeue                                                 
#SBATCH --open-mode=append                                        

# Load modules.
if [ $# -eq 0 ]; then
    module load 2024
    module load OpenMPI/5.0.3-GCC-13.3.0
fi

# Store name of this job script to be used for result naming.
job=$SLURM_JOB_NAME
echo "Job: $job"

# Determine whether the script is executed locally or on Lisa.
if [ $# -eq 0 ]; then
    echo "~ Running on Snellius.."
    WDIR="${HOME}/ComplexityContesters/paper_PF_hardness"
else
    echo "~ Running locally.."
    WDIR="."
fi

# Echo place where the results will finally be stored.
echo "Setting WDIR=$WDIR"

# Compile code to run.
EXP_PATH="${WDIR}/code/experiments"
BIN_NAME="dfs_bnb_mpi"
echo "Going to run $EXP_PATH/$BIN_NAME"

# Get length and set dimension to fold in.
LEN="${job: -2}"
DIM=2

# Set paths for results and data. Then set the dimension to fold.
FIN_RES_PATH="${WDIR}/jobs/${job}_results"
DATASET_PATH="${WDIR}/data/vanEck_random/"
echo "FIN_RES_PATH: $FIN_RES_PATH"
echo "DATASET_PATH: $DATASET_PATH"
mkdir -p "$FIN_RES_PATH"

# Enable checkpointing for prospr's depth_first_bnb algorithm.
export PROSPR_CACHE_DIR=$HOME/prospr_random.cache

# Add the SLURM job id in the results folder.
echo "$SLURM_JOB_ID" >> "$FIN_RES_PATH/slurm_job_ids"
srun "${EXP_PATH}/${BIN_NAME}" "${job}" \
    "${FIN_RES_PATH}" "${DATASET_PATH}" "NULL" "$DIM" "$LEN"

# Requeue job from checkpoint after receiving SIGINT.
if [[ $? -eq 130 ]]; then     
    echo "Received SIGINT, requeueing job.."
    scontrol requeue "$SLURM_JOB_ID"
fi

