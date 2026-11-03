#!/bin/bash
#SBATCH --job-name=md_charmm_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=96:00:00
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

# Load required GROMACS module
module purge
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0

TPR_FILE="step7_production.tpr"
CPT_FILE="step7_production.cpt"
LOG_FILE="step7_production.log"

# Dynamically extract target steps directly from the TPR file
TARGET_STEP=$(gmx dump -s ${TPR_FILE} 2>&1 | awk '/nsteps/ {print $2}' | tr -dc '0-9')
TARGET_STEP=${TARGET_STEP:-500000000}

# Run GROMACS mdrun
gmx_mpi mdrun \
    -s ${TPR_FILE} \
    -cpi ${CPT_FILE} \
    -deffnm step7_production \
    -noappend \
    -ntomp 8 \
    -pin on \
    -nb gpu \
    -bonded gpu \
    -pme gpu \
    -update cpu \
    -cpt 15 \
    -maxh 95.5

# Extract last written step from log file
if [ -f "$LOG_FILE" ]; then
    LAST_STEP=$(grep -E "Step[[:space:]]+Time" -A 1 "$LOG_FILE" | tail -n 1 | awk '{print $1}' | tr -dc '0-9')
fi

# Fallback to checkpoint file if log parsing fails
if [ -z "$LAST_STEP" ]; then
    LAST_STEP=$(gmx check -f ${CPT_FILE} 2>&1 | awk '/^Step/ {print $2}' | tr -dc '0-9')
fi

LAST_STEP=${LAST_STEP:-0}

echo "=========================================="
echo "Current Step: ${LAST_STEP} / ${TARGET_STEP}"
echo "=========================================="

# Resubmission logic
if [ "${LAST_STEP}" -lt "${TARGET_STEP}" ]; then
    echo "Simulation incomplete. Submitting continuation job..."
    sbatch "$0"
else
    echo "Target step ${TARGET_STEP} reached successfully! Stopping pipeline."
fi