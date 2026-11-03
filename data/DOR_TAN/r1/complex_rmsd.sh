#!/bin/bash
#SBATCH --job-name=tan_rmsd
#SBATCH --output=rmsd_%j.out
#SBATCH --error=rmsd_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# ==============================================================================
# SLURM Script: Protein-Ligand (TAN) Combined RMSD Calculation
# ==============================================================================

# Load GROMACS module matching your HPC environment
module purge
module load GROMACS/2025.4-foss-2025b-CUDA-12.9.1

# Input variables
TPR="tan_md_1000ns.tpr"
TRAJ_CLEAN="traj_final_clean_tan_r1.xtc"

echo "=================================================="
echo " Starting RMSD Analysis on $(hostname)"
echo " Date: $(date)"
echo "=================================================="

# Step 1: Create custom index group merging 'Protein' and 'TAN'
echo "--> Creating custom index file for Protein-TAN complex..."
echo -e "\"Protein\" | \"TAN\"\nq" | gmx make_ndx -f $TPR -o index_complex.ndx

# Step 2: Calculate RMSD using the created 'Protein_TAN' group
echo "--> Calculating Protein-TAN Combined RMSD..."
echo -e "Protein_TAN\nProtein_TAN" | gmx rms -s $TPR -f $TRAJ_CLEAN -n index_complex.ndx -o rmsd_complex.xvg -tu ns

echo "=================================================="
echo " Job Finished Successfully at $(date)"
echo " Output written to: rmsd_complex.xvg"
echo "=================================================="