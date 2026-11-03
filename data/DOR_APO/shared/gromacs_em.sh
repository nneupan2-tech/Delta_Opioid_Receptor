#!/bin/bash
#SBATCH -J c1c2b_bryo_restart
#SBATCH -t 1:59:59
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p gpu
#SBATCH --mail-user=nneupan2@cougarnet.uh.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -A DAS

module add GROMACS/2025.2-foss-2023b-CUDA-12.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

gmx_mpi mdrun \
  -deffnm delta_apo_em