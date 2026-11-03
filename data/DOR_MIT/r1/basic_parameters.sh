#!/bin/bash
# ==============================================================================
# Pipeline: Trajectory Concatenation + PBC Fix + Complex Analysis (Protein + Ligand)
# Usage: bash run_complex_analysis.sh
# ==============================================================================

# Input definitions
TPR="mit_md_1000ns.tpr"               # TPR file
TRAJ_FULL="full_concatenated.xtc"     # Concatenated raw trajectory
TRAJ_CLEAN="traj_final_clean.xtc"     # Processed & cleaned output trajectory

# Index groups (Adjust if your index numbers differ):
# Group 1  = Protein
# Group 13 = Ligand (UNK)
# Group 0  = System
# Group 12 = Protein-Ligand Complex (or create dynamically)

echo "=================================================="
echo " STEP 1: CONCATENATING TRAJECTORY PARTS"
echo "=================================================="

# Read all part files safely in numerical order
PARTS=$(ls mit_md_1000ns*.xtc | grep -v "$TRAJ_FULL" | grep -v "$TRAJ_CLEAN" | sort -V)
echo "Found trajectory parts to concatenate:"
echo "$PARTS"

gmx trjcat -f $PARTS -o $TRAJ_FULL

echo "=================================================="
echo " STEP 2: FIXING PBC & CLUSTERING FULL TRAJECTORY"
echo "=================================================="

# 2a. Make all molecules whole across periodic boundaries
echo "--> Making molecules whole..."
echo "0" | gmx trjconv -s $TPR -f $TRAJ_FULL -o step1_whole.xtc -pbc whole

# 2b. Cluster protein and ligand together (Group 1 = Protein, Group 0 = System)
echo "--> Clustering protein-ligand complex together..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step1_whole.xtc -o step2_clustered.xtc -pbc cluster

# 2c. Center clustered complex in simulation box (Group 1 = Protein, Group 0 = System)
echo "--> Centering protein in box..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step2_clustered.xtc -o step3_centered.xtc -center -pbc mol

# 2d. Rotational and translational fitting (Group 1 = Protein, Group 0 = System)
echo "--> Applying rotational & translational fit..."
echo -e "1\n0" | gmx trjconv -s $TPR -f step3_centered.xtc -o $TRAJ_CLEAN -fit rot+trans

# Remove intermediate trajectories to free up disk space
echo "--> Removing temporary intermediate trajectories..."
rm -f $TRAJ_FULL step1_whole.xtc step2_clustered.xtc step3_centered.xtc

echo "=================================================="
echo " STEP 3: RUNNING ANALYSIS ON COMPLEX TRAJECTORY"
echo "=================================================="

# 3a. RMSD Calculations
echo "3a-1. Calculating Protein RMSD..."
echo -e "1\n1" | gmx rms -s $TPR -f $TRAJ_CLEAN -o rmsd_protein.xvg -tu ns

echo "3a-2. Calculating Ligand (UNK) RMSD (Fitted to Protein)..."
echo -e "1\n13" | gmx rms -s $TPR -f $TRAJ_CLEAN -o rmsd_ligand.xvg -tu ns

echo "3a-3. Calculating Protein-Ligand Complex RMSD..."
echo -e "13\n13" | gmx rms -s $TPR -f $TRAJ_CLEAN -o rmsd_complex.xvg -tu ns

# 3b. Protein RMSF
echo "3b. Calculating Protein RMSF..."
echo "1" | gmx rmsf -s $TPR -f $TRAJ_CLEAN -o rmsf_protein.xvg -res

# 3c. Radius of Gyration (Rg)
echo "3c-1. Calculating Protein Radius of Gyration..."
echo "1" | gmx gyrate -s $TPR -f $TRAJ_CLEAN -o gyrate_protein.xvg

echo "3c-2. Calculating Ligand Radius of Gyration..."
echo "13" | gmx gyrate -s $TPR -f $TRAJ_CLEAN -o gyrate_ligand.xvg

# 3d. Protein <-> Ligand (UNK) Hydrogen Bonds (Group 1 vs Group 13)
echo "3d. Calculating Protein-Ligand (UNK) Hydrogen Bonds..."
echo -e "1\n13" | gmx hbond -s $TPR -f $TRAJ_CLEAN -num hbond_prot_unk.xvg -dist hbond_prot_unk_dist.xvg

echo "=================================================="
echo " ALL PROCESSES COMPLETED SUCCESSFULLY!"
echo " Continuous Clean Trajectory: $TRAJ_CLEAN"
echo "=================================================="