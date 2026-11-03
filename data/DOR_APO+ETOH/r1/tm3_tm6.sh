#!/bin/bash
#SBATCH --job-name=tm3_tm6_dist
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --output=tm3_tm6_%j.out
#SBATCH --error=tm3_tm6_%j.err

# ==============================
# 1. Load GROMACS
# ==============================
module purge
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0

# ==============================
# 2. Go to working directory
# ==============================
cd $SLURM_SUBMIT_DIR

# ==============================
# 3. Calculate TM3–TM6 COM distance
# ==============================
gmx_mpi distance \
    -s step7_production.tpr \
    -f traj_final_clean.xtc \
    -n tm3_tm6.ndx \
    -select 'com of group 24 plus com of group 25' \
    -oav tm3_tm6_dist.xvg \
    -tu ns