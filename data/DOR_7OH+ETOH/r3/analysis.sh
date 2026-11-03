#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

echo "=================================================="
echo " STEP 0: CONCATENATING TRAJECTORY CHUNKS"
echo "=================================================="
gmx trjcat -f step7_production.part0001.xtc \
              step7_production.part0002.xtc \
              step7_production.part0003.xtc \
              step7_production.part0004.xtc \
              step7_production.part0005.xtc \
              step7_production.part0006.xtc \
              step7_production.part0007.xtc \
              step7_production.part0008.xtc \
             -o step7_production_combined.xtc

echo "=================================================="
echo " STEP 1: FIXING PBC & CLUSTERING MULTIMERIC CHAINS"
echo "=================================================="
echo "--> Making molecules whole..."
echo "0" | gmx trjconv -s step7_production.tpr -f step7_production_combined.xtc -o step1_whole.xtc -pbc whole

echo "--> Clustering protein chains together..."
echo "1 0" | gmx trjconv -s step7_production.tpr -f step1_whole.xtc -o step2_clustered.xtc -pbc cluster

echo "--> Centering protein in box..."
echo "1 0" | gmx trjconv -s step7_production.tpr -f step2_clustered.xtc -o step3_centered.xtc -center -pbc mol

echo "--> Applying rotational and translational fit..."
echo "1 0" | gmx trjconv -s step7_production.tpr -f step3_centered.xtc -o traj_final_clean.xtc -fit rot+trans

echo "--> Cleaning up temporary trajectory files..."
rm -f step1_whole.xtc step2_clustered.xtc step3_centered.xtc

echo "=================================================="
echo " STEP 2: RUNNING ANALYSIS ON CLEANED TRAJECTORY"
echo "=================================================="
echo "2a. Calculating Protein RMSD..."
echo "1 1" | gmx rms -s step7_production.tpr -f traj_final_clean.xtc -o rmsd.xvg -tu ns

echo "2b. Calculating Protein RMSF..."
echo "1" | gmx rmsf -s step7_production.tpr -f traj_final_clean.xtc -o rmsf.xvg -res

echo "2c. Calculating Protein Radius of Gyration..."
echo "1" | gmx gyrate -s step7_production.tpr -f traj_final_clean.xtc -o gyrate.xvg

echo "2d. Calculating Protein-Ligand (UNK) Hydrogen Bonds..."
echo "1 13" | gmx hbond -s step7_production.tpr -f traj_final_clean.xtc -num hbond_prot_unk.xvg -dist hbond_prot_unk_dist.xvg

echo "2e. Calculating Protein-Ethanol (ETOH) Hydrogen Bonds..."
echo "1 14" | gmx hbond -s step7_production.tpr -f traj_final_clean.xtc -num hbond_prot_etoh.xvg -dist hbond_prot_etoh_dist.xvg

echo "2f. Calculating Ligand (UNK) RMSD..."
# Pipe 1 (Protein) for fitting alignment, Pipe 13 (UNK) for RMSD calculation
echo "1 13" | gmx rms -s step7_production.tpr -f traj_final_clean.xtc -o rmsd_unk.xvg -tu ns

echo "=================================================="
echo " ALL PROCESSES COMPLETED SUCCESSFULLY!"
echo "=================================================="