#!/bin/bash
#SBATCH --job-name=ligand_RMSD
#SBATCH --output=ligand_RMSD_%j.out
#SBATCH --error=ligand_RMSD_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00

module load GROMACS/2025.4-foss-2025b-CUDA-12.9.1

# Files
TPR="step7_production.tpr"
XTC="traj_final_clean_7oh_r3.xtc"
NDX="netra_index.ndx"

# Run ligand RMSD
# Group 4  = Backbone (fit)
# Group 13 = UNK ligand (RMSD)

printf "4\n13\n" | gmx rms \
    -s "$TPR" \
    -f "$XTC" \
    -n "$NDX" \
    -o ligand_RMSD_nm.xvg \
    -tu ns \
    -fit rot+trans

# Convert RMSD from nm to Angstrom
awk '
/^[#@]/ {print; next}
{printf "%.3f %.6f\n", $1, $2*10}
' ligand_RMSD_nm.xvg > ligand_RMSD_fitted.xvg

echo "=========================================="
echo "Ligand RMSD calculation completed"
echo "Fit group     : 4 (Backbone)"
echo "Ligand group  : 13 (UNK)"
echo "Output        : ligand_RMSD_A.xvg"
echo "Units         : Angstrom"
echo "=========================================="