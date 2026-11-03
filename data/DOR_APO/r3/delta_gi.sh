#!/bin/bash
#SBATCH --job-name=gmx_alpha5_rmsd
#SBATCH --output=rmsd_%j.out
#SBATCH --error=rmsd_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# Load GROMACS module (adjust module name to match your cluster environment)
module purge
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0

# Print job information for logging
echo "Starting RMSD calculation at $(date) on host $(hostname)"

# Run relative RMSD calculation
echo "20 22" | gmx rms -s 1000ns.pdb -f traj_final_clean_apo_r3.xtc -n delta_gi_index.ndx -o alpha5_relative_rmsd_fixed.xvg

# Verify completion
if [ $? -eq 0 ]; then
    echo "RMSD calculation completed successfully at $(date)"
else
    echo "RMSD calculation failed!"
    exit 1
fi