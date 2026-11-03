import numpy as np
import pandas as pd
import MDAnalysis as mda

# File definitions
INACTIVE_PDB = "8y45.pdb"
ACTIVE_PDB = "1000ns_tan_r3.pdb"
TRAJECTORY = "traj_final_clean_tan_r3.xtc"
OUTPUT_CSV = "D253_displacement_8y45_over_time.csv"

# 1. Load Universes
ref_u = mda.Universe(ACTIVE_PDB)
mobile_u = mda.Universe(ACTIVE_PDB, TRAJECTORY)

# 2. Extract matching C-alpha atoms for alignment (1-to-1 matching)
ref_ca_all = ref_u.select_atoms("protein and name CA")
mobile_ca_all = mobile_u.select_atoms("protein and name CA")

# Find shared residue IDs
ref_resids = set(ref_ca_all.resids)
mobile_resids = set(mobile_ca_all.resids)
common_resids = sorted(list(ref_resids.intersection(mobile_resids)))

# Select matching atoms from the first chain/monomer
ref_align_atoms = []
mobile_align_atoms = []

for r in common_resids:
    r_atoms = ref_u.select_atoms(f"protein and resid {r} and name CA")
    m_atoms = mobile_u.select_atoms(f"protein and resid {r} and name CA")
    if len(r_atoms) > 0 and len(m_atoms) > 0:
        ref_align_atoms.append(r_atoms[0])
        mobile_align_atoms.append(m_atoms[0])

ref_align_group = mda.AtomGroup(ref_align_atoms)
mobile_align_group = mda.AtomGroup(mobile_align_atoms)

print(
    f"Successfully paired {len(ref_align_group)} matching C-alpha atoms for alignment."
)

# Target D253 atoms
ref_d253 = ref_u.select_atoms("protein and resid 253 and name CA")[0]
mobile_d253 = mobile_u.select_atoms("protein and resid 253 and name CA")[0]

# Reference coordinates (centered at origin)
ref_coords = ref_align_group.positions
ref_center = ref_coords.mean(axis=0)
ref_coords_centered = ref_coords - ref_center
ref_d253_centered = ref_d253.position - ref_center


# Kabsch Alignment Algorithm
def align_coordinates(mobile_coords, ref_coords_centered):
    """Calculates optimal rotation matrix via SVD to align mobile to reference."""
    mob_center = mobile_coords.mean(axis=0)
    mob_centered = mobile_coords - mob_center

    # Covariance matrix
    H = mob_centered.T @ ref_coords_centered
    U, S, Vt = np.linalg.svd(H)
    R = U @ Vt

    # Ensure right-handed coordinate system
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = U @ Vt

    return mob_center, R


# 3. Time-series loop
data = []
for ts in mobile_u.trajectory:
    mob_coords = mobile_align_group.positions

    # Find optimal rotation and centroid
    mob_center, R = align_coordinates(mob_coords, ref_coords_centered)

    # Transform D253 position to aligned reference frame
    mob_d253_pos = mobile_d253.position
    aligned_d253_pos = (mob_d253_pos - mob_center) @ R

    # Measure distance relative to inactive static reference position (Angstroms)
    displacement = np.linalg.norm(aligned_d253_pos - ref_d253_centered)

    data.append(
        {
            "frame": ts.frame,
            "time_ns": round(ts.time / 1000.0, 3),
            "displacement_A": round(displacement, 3),
        }
    )

# Save results
df = pd.DataFrame(data)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Success! Saved frame-by-frame displacement to {OUTPUT_CSV}")