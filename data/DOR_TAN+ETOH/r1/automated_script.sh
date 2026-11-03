#!/bin/bash
#SBATCH -J md_charmm_pipeline
#SBATCH -t 95:59:59
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH -p gpu
#SBATCH --mail-user=nneupan2@cougarnet.uh.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -A DAS
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

set -e

# Load module environment
module purge
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Pipeline Configuration
MAXWARN=1
MAXH=95.5
CPT_INTERVAL=15

# Input / Output Configurations
GRO="system_ethanol.gro"
TOP="topol.top"
NDX="netra_index.ndx"
PROD="step7_production"
STATE_FILE="pipeline_state.txt"
LOG="gromacs_run.log"

exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "GROMACS Pipeline — Job $SLURM_JOB_ID started: $(date)"
echo "Node: $SLURMD_NODENAME"
echo "Resources: 1 node | $SLURM_CPUS_PER_TASK CPU cores | 32G Mem | 1 GPU"
echo "============================================================"

check_exit() {
    if [ $1 -ne 0 ]; then
        echo "[ERROR] '$2' failed at $(date). Check $LOG for details."
        exit 1
    fi
    echo "[OK] '$2' completed at $(date)"
}

mark_done() {
    grep -qx "$1" "$STATE_FILE" 2>/dev/null || echo "$1" >> "$STATE_FILE"
}

is_done() {
    grep -qx "$1" "$STATE_FILE" 2>/dev/null
}

# ---------------------------------------------------------
# STEP 6.0: Energy Minimization
# ---------------------------------------------------------
if ! is_done "step6.0"; then
    echo ""
    echo ">>> STEP 6.0: Energy Minimization"
    gmx_mpi grompp -f step6.0_minimization.mdp \
                   -o step6.0_minimization.tpr \
                   -c "$GRO" -r "$GRO" -p "$TOP" -n "$NDX" \
                   -maxwarn $MAXWARN
    check_exit $? "grompp step6.0"

    gmx_mpi mdrun -v -deffnm step6.0_minimization
    check_exit $? "mdrun step6.0 minimization"
    mark_done "step6.0"
else
    echo "[SKIP] step6.0 already completed."
fi

# ---------------------------------------------------------
# STEP 6.1 to 6.6: Equilibration Stages
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
        check_exit $? "grompp step6.$STEP"

        gmx_mpi mdrun -v -deffnm "$DEFFNM"
        check_exit $? "mdrun step6.$STEP equilibration"
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

if is_done "step7_production"; then
    echo "[DONE] Production MD marked complete in $STATE_FILE. Exiting."
    exit 0
fi

if [ ! -f "${PROD}.tpr" ]; then
    LAST_EQ_GRO="step6.6_equilibration.gro"
    gmx_mpi grompp -f step7_production.mdp \
                   -o "${PROD}.tpr" \
                   -c "$LAST_EQ_GRO" -r "$GRO" -p "$TOP" -n "$NDX" \
                   -maxwarn $MAXWARN
    check_exit $? "grompp step7 production"
else
    echo "[INFO] ${PROD}.tpr already exists — skipping grompp."
fi

# Run or Resume Production MD
set +e
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
set -e

# ---------------------------------------------------------
# VERIFICATION & AUTOMATED RESUBMISSION LOGIC
# ---------------------------------------------------------
TARGET_STEPS=$(grep -i "^nsteps" step7_production.mdp | awk '{print $3}' | tr -d ' \r')
LAST_STEP=$(gmx_mpi check -f "${PROD}.cpt" 2>&1 | grep -i "Step" | awk '{print $2}')

echo "[CHECK] Current Checkpoint Step: $LAST_STEP / Target Step: $TARGET_STEPS"

if [ -n "$LAST_STEP" ] && [ "$LAST_STEP" -ge "$TARGET_STEPS" ]; then
    echo "============================================================"
    echo "[SUCCESS] Target steps ($TARGET_STEPS) reached!"
    echo "SIMULATION FULLY COMPLETE AT: $(date)"
    echo "============================================================"
    mark_done "step7_production"
    exit 0
else
    echo ""
    echo "[RESUBMIT] Simulation incomplete ($LAST_STEP / $TARGET_STEPS steps)."
    echo "[RESUBMIT] Submitting continuation job to Slurm..."
    
    NEW_JOB_ID=$(sbatch --parsable "$0")
    
    echo "[INFO] Submitted continuation Job ID: $NEW_JOB_ID at $(date)"
    echo "============================================================"
    exit 0
fi