#!/bin/bash
#SBATCH --job-name=gmx_hbond
#SBATCH --output=hbond_%j.out
#SBATCH --error=hbond_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=standard     # Adjust to match your partition name

# 1. Purge modules and load missing zlib dependency first
module purge
module load zlib/.1.3.1
module load GROMACS/2025.4-foss-2025b-CUDA-12.9.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 2. Use gmx hbond-legacy for -hbm and -hbn support
echo "1 13" | gmx hbond-legacy \
  -s step7_production.tpr \
  -f traj_final_clean_7oh_r3.xtc \
  -hbm hbmap.xpm \
  -hbn hbond.ndx