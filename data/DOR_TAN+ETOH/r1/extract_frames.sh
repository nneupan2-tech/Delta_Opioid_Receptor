#!/bin/bash
#SBATCH --job-name=extract_frames
#SBATCH --output=extract_frames_%j.out
#SBATCH --error=extract_frames_%j.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# ==================================================
# USER SETTINGS: Change filenames as needed
# ==================================================
TPR="step7_production.tpr"
XTC="traj_final_clean.xtc"

module load GROMACS/2025.4-foss-2025b-CUDA-12.9.1
# ==================================================
# Check files and GROMACS
# ==================================================
set -euo pipefail

command -v gmx >/dev/null 2>&1 || {
    echo "ERROR: gmx not found. Load the appropriate GROMACS module."
    exit 1
}

[[ -f "$TPR" ]] || { echo "ERROR: Missing $TPR"; exit 1; }
[[ -f "$XTC" ]] || { echo "ERROR: Missing $XTC"; exit 1; }

OUTDIR="extracted_frames"
mkdir -p "$OUTDIR"

echo "Starting extraction: $(date)"
echo "TPR: $TPR"
echo "XTC: $XTC"

# Extract frames; GROMACS time units are ps
for NS in 0 200 400 600 800 1000
do
    PS=$((NS * 1000))
    OUTFILE="${OUTDIR}/frame_${NS}ns.gro"

    echo "Extracting ${NS} ns (${PS} ps)..."

    printf "0\n" | gmx trjconv \
        -s "$TPR" \
        -f "$XTC" \
        -o "$OUTFILE" \
        -dump "$PS"

    [[ -s "$OUTFILE" ]] || {
        echo "ERROR: Output missing or empty: $OUTFILE"
        exit 1
    }

    echo "Saved: $OUTFILE"
done

echo "Extraction completed: $(date)"
