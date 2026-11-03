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
# 1. Load Modules
# ==============================
module purge
module load zlib/.1.2.13-GCCcore-13.2.0
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0

# ==============================
# 2. Environment Fixes
# ==============================
# Bypass PSM2/Omni-Path fabric errors on single-node jobs
export OMPI_MCA_mtl=^psm2
export OMPI_MCA_pml=ob1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ==============================
# 3. Go to working directory
# ==============================
cd $SLURM_SUBMIT_DIR

# ==============================
# 4. Calculate TM3–TM6 COM distance
# ==============================
gmx_mpi distance \
    -s apo_md_1000ns.tpr \
    -f traj_final_clean_apo_r1.xtc \
    -n tm3_tm6.ndx \
    -select 'com of group 23 plus com of group 24' \
    -oav tm3_tm6_dist.xvg \
    -tu ns