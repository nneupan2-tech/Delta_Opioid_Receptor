#!/bin/bash
#SBATCH --job-name=mmpbsa_run
#SBATCH --nodes=1
#SBATCH --ntasks=8            
#SBATCH --mem=64G              
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --output=job_output_%j.log
#SBATCH --error=job_error_%j.err

# 1. Environment setup
module purge
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0
source activate /project/das/Netra/conda_envs/gmx_MMPBSA

# 2. Run execution with 8 MPI tasks
mpirun -np $SLURM_NTASKS gmx_MMPBSA \
  -O \
  -i mmpbsa.in \
  -cs tan_md_1000ns.tpr \
  -cp topol.top \
  -ct traj_final_clean_tan_r3.xtc \
  -ci netra_index.ndx \
  -cg 1 13 \
  -o FINAL_RESULTS_MMPBSA.dat \
  -eo PER_RESIDUE_ENERGIES.csv \
  -deo FINAL_DECOMP_MMPBSA.dat