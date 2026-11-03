#!/bin/bash
#SBATCH -J md_charmm_pipeline
#SBATCH -t 95:59:59
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p gpu
#SBATCH --mail-user=nneupan2@cougarnet.uh.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -A DAS
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

module add GROMACS/2025.2-foss-2023b-CUDA-12.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

MAXH=95.5
CPT_INTERVAL=15
PROD="step7_production"
STATE_FILE="pipeline_state.txt"
LOG="gromacs_run.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "GROMACS Continuation — Job $SLURM_JOB_ID started: $(date)"
echo "Node: $SLURMD_NODENAME"
echo "Resources: 1 node | $SLURM_CPUS_PER_TASK CPU cores | 1 GPU"
echo "============================================================"

# Mark earlier steps complete so it skips minimization & equilibration
echo "step6.0" >> "$STATE_FILE"
for s in 1 2 3 4 5 6; do echo "step6.$s" >> "$STATE_FILE"; done

# Check if production is already done
if [ -f "${PROD}.gro" ]; then
    echo "[SUCCESS] Production MD is already complete! (${PROD}.gro found)"
    exit 0
fi

if [ ! -f "${PROD}.cpt" ]; then
    echo "[ERROR] ${PROD}.cpt not found! Cannot resume from checkpoint."
    exit 1
fi

echo "[INFO] Resuming Production MD from checkpoint (${PROD}.cpt)..."

# Temporarily disable 'exit on error' so SLURM walltime timeout doesn't crash the script before resubmission
set +e

gmx_mpi mdrun \
    -s "${PROD}.tpr" \
    -cpi "${PROD}.cpt" \
    -deffnm "$PROD" \
    -noappend \
    -ntomp $SLURM_CPUS_PER_TASK \
    -pin on \
    -nb gpu \
    -bonded gpu \
    -pme gpu \
    -update cpu \
    -cpt $CPT_INTERVAL \
    -maxh $MAXH

MDRUN_EXIT_CODE=$?
set -e

# Check if simulation finished completely during this run
if [ -f "${PROD}.gro" ]; then
    echo "============================================================"
    echo "SIMULATION FULLY COMPLETED: $(date)"
    echo "============================================================"
    exit 0
fi

# If .gro is still missing, walltime was reached — submit the next continuation job automatically
echo ""
echo "[RESUBMIT] Walltime reached ($MAXH hours). Resubmitting continuation job to Sabine..."
NEW_JOB_ID=$(sbatch --parsable "$0")
echo "[INFO] Next SLURM Job ID: $NEW_JOB_ID submitted at $(date)"
echo "[INFO] Track status with: squeue -u $USER"

exit 0