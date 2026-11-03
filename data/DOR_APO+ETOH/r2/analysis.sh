#!/bin/bash
# ==============================================================================
# Master Pipeline: Trajectory Concatenation, PBC Cleanup & Full Structural Analysis
# Usage: bash fix_and_analyze.sh
# ==============================================================================

# Input trajectory & reference files
TPR="step7_production.tpr"
TRAJ_COMBINED="step7_production_combined.xtc"
TRAJ_CLEAN="traj_final_clean.xtc"
INDEX_FILE="hbond.ndx"

echo "=================================================="
echo " STEP 0: CONCATENATING MULTIPLE TRAJECTORY CHUNKS"
echo "=================================================="

# 0a. Combine part0001, part0002, and part0003 into one trajectory
echo "--> Combining trajectory parts (part0001, part0002, part0003)..."
gmx trjcat -f step7_production.part0001.xtc \
              step7_production.part0002.xtc \
              step7_production.part0003.xtc \
              step7_production.part0004.xtc \
              step7_production.part0005.xtc \
              step7_production.part0006.xtc \
              step7_production.part0007.xtc \
              step7_production.part0008.xtc \
           -o $TRAJ_COMBINED

echo "=================================================="
echo " STEP 1: FIXING PBC & CLUSTERING MULTIMERIC CHAINS"
echo "=================================================="

# 1a. Make molecules whole across periodic boundaries
echo "--> Making molecules whole..."
echo "0" | gmx trjconv -s $TPR -f $TRAJ_COMBINED -o step1_whole.xtc -pbc whole

# 1b. Cluster all protein chains together (fixes multi-chain boundary splits)
# Group 1 = Protein (for clustering), Group 0 = System (output)
echo "--> Clustering protein chains together..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step1_whole.xtc -o step2_clustered.xtc -pbc cluster

# 1c. Center the clustered protein in the simulation box
# Group 1 = Protein (centering), Group 0 = System (output)
echo "--> Centering protein in box..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step2_clustered.xtc -o step3_centered.xtc -center -pbc mol

# 1d. Rotational and translational fitting on the clean protein
# Group 1 = Protein (fitting), Group 0 = System (output)
echo "--> Applying rotational and translational fit..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step3_centered.xtc -o $TRAJ_CLEAN -fit rot+trans

# Remove intermediate files to save disk space
echo "--> Cleaning up intermediate trajectories..."
rm -f step1_whole.xtc step2_clustered.xtc step3_centered.xtc $TRAJ_COMBINED

echo "=================================================="
echo " STEP 2: RUNNING ANALYSIS ON CLEANED TRAJECTORY"
echo "=================================================="

# 2a. Protein RMSD (Group 1 vs Group 1)
echo "2a. Calculating Protein RMSD..."
echo -e "1\n1" | gmx rms -s $TPR -f $TRAJ_CLEAN -o rmsd.xvg -tu ns

# 2b. Protein RMSF per residue (Group 1)
echo "2b. Calculating Protein RMSF..."
echo "1" | gmx rmsf -s $TPR -f $TRAJ_CLEAN -o rmsf.xvg -res

# 2c. Radius of Gyration (Group 1)
echo "2c. Calculating Protein Radius of Gyration..."
echo "1" | gmx gyrate -s $TPR -f $TRAJ_CLEAN -o gyrate.xvg

# 2d. Protein <-> Ethanol (ETOH) Hydrogen Bonds
echo "2d. Calculating Protein-Ethanol (ETOH) Hydrogen Bonds..."
if [ -f "$INDEX_FILE" ]; then
    echo -e "Protein\nETOH" | gmx hbond -s $TPR -f $TRAJ_CLEAN -n $INDEX_FILE -num hbond_prot_etoh.xvg -dist hbond_prot_etoh_dist.xvg
else
    echo -e "1\n19" | gmx hbond -s $TPR -f $TRAJ_CLEAN -num hbond_prot_etoh.xvg -dist hbond_prot_etoh_dist.xvg
fi

echo "=================================================="
echo " ALL PROCESSES COMPLETED SUCCESSFULLY!"
echo "--------------------------------------------------"
echo " Cleaned Trajectory: $TRAJ_CLEAN"
echo " Output Analysis Files:"
echo "   - RMSD:                  rmsd.xvg"
echo "   - RMSF:                  rmsf.xvg"
echo "   - Radius of Gyration:    gyrate.xvg"
echo "   - Protein-Ethanol HBonds:hbond_prot_etoh.xvg"
echo "=================================================="