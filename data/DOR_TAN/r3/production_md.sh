#!/bin/bash
#SBATCH -J tan_gmx_production
#SBATCH --mail-user=nneupan2@cougarnet.uh.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -N 1
#SBATCH -t 3-23:50:00
#SBATCH --ntasks=8
#SBATCH --gres=gpu:1 
#SBATCH --mem-per-cpu=16GB
#SBATCH -o slurm_%j.log
#SBATCH -e slurm_%j.err

set -euo pipefail

# Load GROMACS environment
module purge
module add GROMACS/2025.2-foss-2023b-CUDA-12.4.0

# ====================================================================
# CONFIGURATION: CHANGE THIS NAME PER FOLDER TO MATCH YOUR .TPR FILE
# ====================================================================
DEFFNM="tan_md_1000ns"
DONE_FILE=".md_done"

echo "Job ${SLURM_JOB_ID} started at $(date)"
echo "Running in $(pwd)"

# 1. Safeguards
if [ ! -f "${DEFFNM}.tpr" ]; then
    echo "ERROR: ${DEFFNM}.tpr not found!"
    exit 1
fi

if [ -f "$DONE_FILE" ]; then
    echo "Simulation already marked complete. No action taken."
    exit 0
fi

# 2. Run GROMACS (Auto-Detect Fresh vs Resume)
if [ -f "${DEFFNM}.cpt" ]; then
    echo "[INFO] Checkpoint found (${DEFFNM}.cpt). Resuming simulation..."
    gmx_mpi mdrun \
        -deffnm "${DEFFNM}" \
        -cpi "${DEFFNM}.cpt" \
        -noappend \
        -maxh 95.5
else
    echo "[INFO] No checkpoint found. Starting simulation fresh from TPR..."
    gmx_mpi mdrun \
        -deffnm "${DEFFNM}" \
        -noappend \
        -maxh 95.5
fi

# 3. Dynamic Continuation Check
if [ -f "${DEFFNM}.gro" ]; then
    echo "Simulation reached the absolute final step successfully! (${DEFFNM}.gro found)"
    touch "$DONE_FILE"
else
    echo "Simulation timed out (-maxh reached). Auto-submitting next chunk..."
    NEXT_JOB=$(sbatch --parsable "$0")
    echo "Next job submitted successfully with ID: $NEXT_JOB"
fi

echo "Job ${SLURM_JOB_ID} ended at $(date)"