#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# This script folds the proteins of length 25 for the hratio dataset.
# Pass an argument, doesn't matter what, to trigger local execution.
# The default is Snellius execution.
# With 128 tasks and 1 node, this script takes around ?? minutes to run.
# ------------------------------------------------------------------------------
#SBATCH --job=hratio_l25
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
TMP_RES_PATH="${TMPDIR}/${job}_results"
FIN_RES_PATH="${WDIR}/jobs/${job}_results"
DATASET_PATH="${WDIR}/data/vanEck_hratio/$LEN"
echo "TMP_RES_PATH: $TMP_RES_PATH"
echo "FIN_RES_PATH: $FIN_RES_PATH"
echo "DATASET_PATH: $DATASET_PATH"
mkdir -p "$FIN_RES_PATH"
mkdir -p "$TMP_RES_PATH"

# Add the SLURM job id in the results folder.
echo "$SLURM_JOB_ID" >> "$FIN_RES_PATH/slurm_job_ids"

# Copy existing results if there are any.
if [[ -d "$FIN_RES_PATH" ]]; then
    rsync -au "$FIN_RES_PATH/" "$TMP_RES_PATH"
fi

# ChatGPT sed stuff for extracting the H<num> part of the data files.
h_values=$(ls "${DATASET_PATH}" | sed -E 's/.*_(H[0-9]+)\.csv/\1/' | sort -u)

# Loop over all hratio options and solve one by one.
for h in $h_values; do
    # Compute remaining time minus 2 minutes for copying results.                
    time_left=$(squeue -j "$SLURM_JOB_ID" -h -o "%L")                            
    
    # First extract days if present.
    if [[ "$time_left" == *-* ]]; then
        days=${time_left%%-*}
        time_part=${time_left#*-}
    else
        days=0
        time_part=$time_left
    fi

    # Then compute the time.
    IFS=':' read -ra t <<< "$time_part"
    case ${#t[@]} in                                                             
        3) total_seconds=$((10#$days*86400 + 10#${t[0]}*3600 + 10#${t[1]}*60 + 10#${t[2]})) ;;    
        2) total_seconds=$((10#$days*86400 + 10#${t[0]}*60 + 10#${t[1]})) ;;                      
        1) total_seconds=$((10#$days*86400 + 10#${t[0]})) ;;                                      
    esac                                                                         
    
    # Remove 3 minutes for copying over the results and shutting down.
    timeout_seconds=$(( total_seconds > 180 ? total_seconds - 180 : total_seconds ))

    echo "Processing $h.."
    timeout "${timeout_seconds}s" srun "${EXP_PATH}/${BIN_NAME}" "${job}" \
        "${TMP_RES_PATH}" "${DATASET_PATH}" "$h" $DIM $LEN

    # Copy over results.
    rsync -au "${TMP_RES_PATH}/" "${FIN_RES_PATH}"
done
