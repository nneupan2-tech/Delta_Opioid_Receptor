#!/bin/bash
# ==============================================================================
# Pipeline: Trajectory Concatenation + PBC Fix + Analysis for Multi-part Runs
# Usage: bash run_apo_analysis.sh
# ==============================================================================

# Input definitions
TPR="apo_md_1000ns.tpr"               # Adjust filename if your .tpr has a different name
TRAJ_FULL="full_concatenated.xtc"     # Concatenated raw trajectory
TRAJ_CLEAN="traj_final_clean.xtc"     # Processed & cleaned output trajectory

echo "=================================================="
echo " STEP 1: CONCATENATING TRAJECTORY PARTS (1 to 6)"
echo "=================================================="

# gmx trjcat will read all parts, remove duplicates, and concatenate in order
gmx trjcat -f apo_md_1000ns.part*.xtc -o $TRAJ_FULL

echo "=================================================="
echo " STEP 2: FIXING PBC & CLUSTERING FULL TRAJECTORY"
echo "=================================================="

# 2a. Make all molecules whole across periodic boundaries
echo "--> Making molecules whole..."
echo "0" | gmx trjconv -s $TPR -f $TRAJ_FULL -o step1_whole.xtc -pbc whole

# 2b. Cluster protein chains together (Group 1 = Protein, Group 0 = System)
echo "--> Clustering protein chains together..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step1_whole.xtc -o step2_clustered.xtc -pbc cluster

# 2c. Center clustered protein in simulation box (Group 1 = Protein, Group 0 = System)
echo "--> Centering protein in box..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step2_clustered.xtc -o step3_centered.xtc -center -pbc mol

# 2d. Rotational and translational fitting (Group 1 = Protein, Group 0 = System)
echo "--> Applying rotational & translational fit..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step3_centered.xtc -o $TRAJ_CLEAN -fit rot+trans

# Remove intermediate trajectories to free up disk space
echo "--> Removing temporary intermediate trajectories..."
rm -f $TRAJ_FULL step1_whole.xtc step2_clustered.xtc step3_centered.xtc

echo "=================================================="
echo " STEP 3: RUNNING ANALYSIS ON FULL 1000 ns TRAJECTORY"
echo "=================================================="

# 3a. Protein RMSD
echo "3a. Calculating Protein RMSD..."
echo -e "1\n1" | gmx rms -s $TPR -f $TRAJ_CLEAN -o rmsd.xvg -tu ns

# 3b. Protein RMSF
echo "3b. Calculating Protein RMSF..."
echo "1" | gmx rmsf -s $TPR -f $TRAJ_CLEAN -o rmsf.xvg -res

# 3c. Radius of Gyration (Rg)
echo "3c. Calculating Protein Radius of Gyration..."
echo "1" | gmx gyrate -s $TPR -f $TRAJ_CLEAN -o gyrate.xvg

# 3d. Protein <-> Ligand (UNK) Hydrogen Bonds (Group 1 vs Group 13)
# Note: If this is an APO system without UNK/Ligand, you can comment this block out
echo "3d. Calculating Protein-Ligand (UNK) Hydrogen Bonds..."
echo -e "1\n13" | gmx hbond -s $TPR -f $TRAJ_CLEAN -num hbond_prot_unk.xvg -dist hbond_prot_unk_dist.xvg

# 3e. Protein <-> Ethanol (ETOH) Hydrogen Bonds (Group 1 vs Group 20)
# Note: Adjust group index '20' if ETOH index changed in this TPR file
echo "3e. Calculating Protein-Ethanol (ETOH) Hydrogen Bonds..."
echo -e "1\n20" | gmx hbond -s $TPR -f $TRAJ_CLEAN -num hbond_prot_etoh.xvg -dist hbond_prot_etoh_dist.xvg

echo "=================================================="
echo " ALL PROCESSES COMPLETED SUCCESSFULLY!"
echo " Continuous Clean Trajectory: $TRAJ_CLEAN"
echo "=================================================="