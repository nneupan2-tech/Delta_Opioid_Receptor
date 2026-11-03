#!/bin/bash
#SBATCH --job-name=rmsds_fix_r2
#SBATCH --output=rmsds_fix_%j.out
#SBATCH --error=rmsds_fix_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# #SBATCH --partition=normal    # Uncomment and edit if your cluster requires a partition

module load GROMACS/2025.4-foss-2025b-CUDA-12.9.1
cd $SLURM_SUBMIT_DIR

TPR="tan_md_1000ns.tpr"
TRAJ_IN="traj_final_clean_tan_r1.xtc"
TRAJ_OUT="traj_final_fixed_tan_r1.xtc"

echo "=================================================="
echo " Starting Fast PBC Fix: Ligand & Complex RMSD"
echo " Job ID: $SLURM_JOB_ID"
echo "=================================================="

# 1. Create combined index group for Protein + TAN
echo "--> Creating index group for Protein_TAN..."
echo -e "\"Protein\" | \"TAN\"\nq" | gmx make_ndx -f $TPR -o index_complex.ndx

# 2. Fix PBC wrapping on existing trajectory
echo "--> Re-clustering ligand back to protein..."
echo -e "Protein_TAN\n0" | gmx trjconv -s $TPR -f $TRAJ_IN -n index_complex.ndx -o temp1.xtc -pbc cluster
echo -e "Protein_TAN\n0" | gmx trjconv -s $TPR -f temp1.xtc -n index_complex.ndx -o temp2.xtc -center -pbc mol -ur compact
echo -e "1\n0" | gmx trjconv -s $TPR -f temp2.xtc -o $TRAJ_OUT -fit rot+trans

rm -f temp1.xtc temp2.xtc

# 3a. Calculate Ligand RMSD (Fit to Group 1: Protein, Calculate on Group 13: TAN)
echo "--> Calculating Ligand RMSD..."
echo -e "1\n13" | gmx rms -s $TPR -f $TRAJ_OUT -o rmsd_ligand_tan_r2_fixed.xvg -tu ns

# 3b. Calculate Complex RMSD (Fit to Protein_TAN, Calculate on Protein_TAN)
echo "--> Calculating Protein-Ligand Complex RMSD..."
echo -e "Protein_TAN\nProtein_TAN" | gmx rms -s $TPR -f $TRAJ_OUT -n index_complex.ndx -o rmsd_complex_tan_r2_fixed.xvg -tu ns

echo "=================================================="
echo " Job Complete!"
echo " Outputs saved as:"
echo " - rmsd_ligand_tan_r2_fixed.xvg"
echo " - rmsd_complex_tan_r2_fixed.xvg"
echo "=================================================="