#!/bin/bash
#SBATCH -J mit_md_production
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

# 1. Start logging immediately
LOG="gromacs_mit_production_run.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "GROMACS Production MD (MIT) — Job $SLURM_JOB_ID started: $(date)"
echo "Node: $SLURMD_NODENAME"
echo "============================================================"

# 2. Turn on strict error checking
set -e

# 3. Load GROMACS module (bypass CMake version reloads gracefully)
module add GROMACS/2025.2-foss-2023b-CUDA-12.4.0 || true
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Core Variables updated for your new system
PROD="mit_md_1000ns"
CPT_INTERVAL=15
MAXH=95.5

# Safety check: If the final production .gro file already exists, stop everything.
if [ -f "${PROD}.gro" ]; then
    echo "[DONE] Production MD already completed (${PROD}.gro found). Exiting."
    exit 0
fi

# 4. IMMEDIATE DEPENDENCY SUBMISSION
# This submits the next cycle immediately to secure your place in the queue.
# It will stay in 'PD' status and ONLY run if this current Job ID finishes successfully.
NEXT_JOB_ID=$(sbatch --parsable --dependency=afterok:$SLURM_JOB_ID "$0")
echo "[DEPENDENCY] Next job $NEXT_JOB_ID submitted. Standing by in queue..."
echo "============================================================"

# 5. Run GROMACS mdrun
if [ -f "${PROD}.cpt" ]; then
    echo "[INFO] Checkpoint found — Resuming production MD from ${PROD}.cpt"
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
else
    echo "[INFO] No checkpoint found — Starting production MD from scratch using ${PROD}.tpr"
    gmx_mpi mdrun \
        -s "${PROD}.tpr" \
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
fi

# 6. Finalization Check
if [ -f "${PROD}.gro" ]; then
    echo "[SUCCESS] Production MD completed fully!"
    echo "============================================================"
    echo "SIMULATION COMPLETE: $(date)"
    echo "============================================================"
    exit 0
fi

echo "[INFO] Clean exit via walltime limit. Handing over to job $NEXT_JOB_ID."
exit 0