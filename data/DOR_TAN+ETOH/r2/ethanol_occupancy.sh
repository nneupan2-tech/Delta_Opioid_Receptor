#!/bin/bash
#SBATCH --job-name=eth_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --output=eth_anal_%j.log
#SBATCH --error=eth_anal_%j.err

# 1. Environment setup
module purge
module load GROMACS/2025.2-foss-2023b-CUDA-12.4.0

# 2. File definitions
TPR="step7_production.tpr"
TRAJ="traj_final_clean.xtc"
NDX="etoh_index.ndx"

echo "=== 1. Calculating Minimum Distance and Contact Counts ==="
# Analyzes minimum distance and contact numbers between Group 13 (Drug) and ETOH
printf "13\n20\n" | gmx mindist \
  -s step7_production.tpr \
  -f traj_final_clean.xtc \
  -n etoh_index.ndx \
  -od mindist_ethanol_drug.xvg \
  -on numcont_ethanol_drug.xvg \
  -d 0.4

echo "=== 2. Calculating Binding Pocket Ethanol Occupancy ==="
# Counts ETOH residues within 0.4 nm (4 Angstroms) of Group 13 (Drug)
gmx select \
  -s ${TPR} \
  -f ${TRAJ} \
  -n ${NDX} \
  -select 'resname ETOH and within 0.4 of group 13' \
  -os ethanol_occupancy.xvg

echo "Analysis complete."