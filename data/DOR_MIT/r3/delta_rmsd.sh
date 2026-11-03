#!/bin/bash
#SBATCH --job-name=dOR_RMSD
#SBATCH --output=dOR_RMSD_%j.out
#SBATCH --error=dOR_RMSD_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --partition=batch

module load GROMACS/2025.4-foss-2025b-CUDA-12.9.1

TPR_FILE="mit_md_1000ns.tpr"
XTC_FILE="traj_final_clean_mit_r3.xtc"

INDEX_FILE="index_dOR_backbone.ndx"
OUT_XVG="dOR_backbone_RMSD.xvg"

echo "=========================================="
echo "d-opioid receptor RMSD"
echo "Chain A, residues 45-330, backbone"
echo "=========================================="

# Create dOR backbone index
gmx select \
    -s "$TPR_FILE" \
    -select 'chain A and resid 45 to 330 and backbone' \
    -on "$INDEX_FILE"

echo ""
echo "Index file created:"
cat "$INDEX_FILE"

echo ""
echo "Calculating dOR backbone RMSD..."

# IMPORTANT:
# Replace 0 0 below with the actual index-group number
# corresponding to the generated dOR selection.
printf "0\n0\n" | gmx rms \
    -s "$TPR_FILE" \
    -f "$XTC_FILE" \
    -n "$INDEX_FILE" \
    -o "$OUT_XVG" \
    -tu ns

rm -f "$INDEX_FILE"

echo ""
echo "=========================================="
echo "Done."
echo "Output: $OUT_XVG"
echo "=========================================="