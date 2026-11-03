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

set -e

# Load the exact GROMACS MPI module from your environment
module add GROMACS/2025.2-foss-2023b-CUDA-12.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Pipeline Configuration
MAXWARN=1
MAXH=95.5
CPT_INTERVAL=15

# Core Input Files (Updated based on your directory)
GRO="system_mixed.gro"
TOP="topol.top"
NDX="netra_index.ndx"
PROD="step7_production"
STATE_FILE="pipeline_state.txt"
LOG="gromacs_run.log"

exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "GROMACS Pipeline — Job $SLURM_JOB_ID started: $(date)"
echo "Node: $SLURMD_NODENAME"
echo "Resources: 1 node | $SLURM_CPUS_PER_TASK CPU cores | 1 GPU"
echo "============================================================"

check_exit() {
    if [ $? -ne 0 ]; then
        echo "[ERROR] '$1' failed at $(date). Check $LOG for details."
        exit 1
    fi
    echo "[OK] '$1' completed at $(date)"
}

mark_done() {
    grep -qx "$1" "$STATE_FILE" 2>/dev/null || echo "$1" >> "$STATE_FILE"
}

is_done() {
    grep -qx "$1" "$STATE_FILE" 2>/dev/null
}

# ---------------------------------------------------------
# STEP 6.0: Minimization
# ---------------------------------------------------------
if ! is_done "step6.0"; then
    echo ""
    echo ">>> STEP 6.0: Energy Minimization"
    gmx_mpi grompp -f step6.0_minimization.mdp \
               -o step6.0_minimization.tpr \
               -c "$GRO" -r "$GRO" -p "$TOP" -n "$NDX" \
               -maxwarn $MAXWARN
    check_exit "grompp step6.0"

    gmx_mpi mdrun -v -deffnm step6.0_minimization
    check_exit "mdrun step6.0 minimization"
    mark_done "step6.0"
else
    echo "[SKIP] step6.0 already completed."
fi

# ---------------------------------------------------------
# STEP 6.1 to 6.6: Equilibration
# ---------------------------------------------------------
for STEP in 1 2 3 4 5 6; do
    DEFFNM="step6.${STEP}_equilibration"
    MDP="step6.${STEP}_equilibration.mdp"

    if ! is_done "step6.$STEP"; then
        echo ""
        echo ">>> STEP 6.$STEP: Equilibration — $MDP"

        if [ $STEP -eq 1 ]; then
            PREV_GRO="step6.0_minimization.gro"
        else
            PREV_GRO="step6.$((STEP-1))_equilibration.gro"
        fi

        gmx_mpi grompp -f "$MDP" \
                   -o "${DEFFNM}.tpr" \
                   -c "$PREV_GRO" -r "$GRO" -p "$TOP" -n "$NDX" \
                   -maxwarn $MAXWARN
        check_exit "grompp step6.$STEP"

        gmx_mpi mdrun -v -deffnm "$DEFFNM"
        check_exit "mdrun step6.$STEP equilibration"
        mark_done "step6.$STEP"
    else
        echo "[SKIP] step6.$STEP already completed."
    fi
done

# ---------------------------------------------------------
# STEP 7: Production MD
# ---------------------------------------------------------
echo ""
echo ">>> STEP 7: Production MD"

if [ ! -f "${PROD}.tpr" ]; then
    LAST_EQ_GRO="step6.6_equilibration.gro"
    gmx_mpi grompp -f step7_production.mdp \
               -o "${PROD}.tpr" \
               -c "$LAST_EQ_GRO" -r "$GRO" -p "$TOP" -n "$NDX" \
               -maxwarn $MAXWARN
    check_exit "grompp step7 production"
else
    echo "[INFO] ${PROD}.tpr already exists — skipping grompp."
fi

if is_done "step7_production"; then
    echo "[DONE] Production MD already completed. Nothing to do."
    exit 0
fi

if [ -f "${PROD}.cpt" ]; then
    echo "[INFO] Checkpoint found — resuming production MD."
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
    echo "[INFO] No checkpoint found — starting production MD from scratch."
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

if [ -f "${PROD}.gro" ]; then
    echo "[SUCCESS] Production MD completed fully!"
    mark_done "step7_production"
    echo "============================================================"
    echo "SIMULATION COMPLETE: $(date)"
    echo "============================================================"
    exit 0
fi

# ---------------------------------------------------------
# AUTO-RESUBMISSION LOGIC
# ---------------------------------------------------------
echo ""
echo "[RESUBMIT] Wall time reached — resubmitting job to Sabine queue..."
NEW_JOB_ID=$(sbatch --parsable "$0")
echo "[INFO] New SLURM Job ID: $NEW_JOB_ID submitted at $(date)"
echo "[INFO] Monitor with: squeue -u $USER"

exit 0