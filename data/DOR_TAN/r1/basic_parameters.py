# -*- coding: utf-8 -*-
import os
import matplotlib

# Force non-interactive backend for HPC cluster execution
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np


def read_xvg(filename):
    """Parses GROMACS .xvg files safely ignoring header lines."""
    data = []
    if not os.path.exists(filename):
        print(f'Warning: File {filename} not found!')
        return None
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                cols = line.strip().split()
                if len(cols) >= 2:
                    data.append([float(cols[0]), float(cols[1])])
    return np.array(data)


# Match exact file names for the 'tan' ligand system
files = {
    'rmsd_protein': ['rmsd_protein_tan_r1.xvg', 'rmsd_protein_tan_r2.xvg', 'rmsd_protein_tan_r3.xvg'],
    'rmsd_ligand':  ['rmsd_ligand_tan_r1.xvg', 'rmsd_ligand_tan_r2.xvg', 'rmsd_ligand_tan_r3.xvg'],
    'gyrate':       ['gyrate_tan_r1.xvg', 'gyrate_tan_r2.xvg', 'gyrate_tan_r3.xvg'],
    'hbond':        ['hbond_prot_unk_tan_r1.xvg', 'hbond_prot_unk_tan_r2.xvg', 'hbond_prot_unk_tan_r3.xvg'],
    'rmsf':         ['rmsf_tan_r1.xvg', 'rmsf_tan_r2.xvg', 'rmsf_tan_r3.xvg']
}

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
labels = ['Replicate 1', 'Replicate 2', 'Replicate 3']

# Figure Setup: Top 2 rows for system metrics, Bottom 2 rows for 4 Chain RMSFs
fig = plt.figure(figsize=(16, 18), dpi=300)
gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 1], hspace=0.35, wspace=0.25)

# Axis definitions
ax_rmsd_prot = fig.add_subplot(gs[0, 0])
ax_rmsd_lig  = fig.add_subplot(gs[0, 1])
ax_gyrate    = fig.add_subplot(gs[1, 0])
ax_hbond     = fig.add_subplot(gs[1, 1])

# Subplots for 4 individual protein chains
ax_chains = [
    fig.add_subplot(gs[2, 0]),
    fig.add_subplot(gs[2, 1]),
    fig.add_subplot(gs[3, 0]),
    fig.add_subplot(gs[3, 1]),
]

# ==========================================
# 1. PROTEIN RMSD
# ==========================================
for i, f_name in enumerate(files['rmsd_protein']):
    data = read_xvg(f_name)
    if data is not None:
        ax_rmsd_prot.plot(data[:, 0], data[:, 1], color=colors[i], label=labels[i], lw=1, alpha=0.85)

ax_rmsd_prot.set_title('A. Protein Backbone RMSD', fontsize=12, fontweight='bold', loc='left')
ax_rmsd_prot.set_xlabel('Time (ns)', fontsize=10)
ax_rmsd_prot.set_ylabel('RMSD (nm)', fontsize=10)
ax_rmsd_prot.legend(loc='upper left', frameon=True, fontsize=9)
ax_rmsd_prot.grid(True, linestyle='--', alpha=0.5)

# ==========================================
# 2. LIGAND RMSD
# ==========================================
for i, f_name in enumerate(files['rmsd_ligand']):
    data = read_xvg(f_name)
    if data is not None:
        ax_rmsd_lig.plot(data[:, 0], data[:, 1], color=colors[i], label=labels[i], lw=1, alpha=0.85)

ax_rmsd_lig.set_title('B. Ligand RMSD', fontsize=12, fontweight='bold', loc='left')
ax_rmsd_lig.set_xlabel('Time (ns)', fontsize=10)
ax_rmsd_lig.set_ylabel('RMSD (nm)', fontsize=10)
ax_rmsd_lig.legend(loc='upper left', frameon=True, fontsize=9)
ax_rmsd_lig.grid(True, linestyle='--', alpha=0.5)

# ==========================================
# 3. RADIUS OF GYRATION (Rg)
# ==========================================
for i, f_name in enumerate(files['gyrate']):
    data = read_xvg(f_name)
    if data is not None:
        time_ns = data[:, 0] / 1000.0 if data[:, 0].max() > 5000 else data[:, 0]
        ax_gyrate.plot(time_ns, data[:, 1], color=colors[i], label=labels[i], lw=1, alpha=0.85)

ax_gyrate.set_title('C. Radius of Gyration ($R_g$)', fontsize=12, fontweight='bold', loc='left')
ax_gyrate.set_xlabel('Time (ns)', fontsize=10)
ax_gyrate.set_ylabel('$R_g$ (nm)', fontsize=10)
ax_gyrate.legend(loc='upper right', frameon=True, fontsize=9)
ax_gyrate.grid(True, linestyle='--', alpha=0.5)

# ==========================================
# 4. HYDROGEN BONDS (Protein-Ligand)
# ==========================================
for i, f_name in enumerate(files['hbond']):
    data = read_xvg(f_name)
    if data is not None:
        time_ns = data[:, 0] / 1000.0 if data[:, 0].max() > 5000 else data[:, 0]
        ax_hbond.plot(time_ns, data[:, 1], color=colors[i], label=labels[i], lw=0.8, alpha=0.7)

ax_hbond.set_title('D. Protein-Ligand Hydrogen Bonds', fontsize=12, fontweight='bold', loc='left')
ax_hbond.set_xlabel('Time (ns)', fontsize=10)
ax_hbond.set_ylabel('Number of H-Bonds', fontsize=10)
ax_hbond.legend(loc='upper right', frameon=True, fontsize=9)
ax_hbond.grid(True, linestyle='--', alpha=0.5)

# ==========================================
# 5. SEPARATED RMSF PER PROTEIN CHAIN
# ==========================================
chain_names = ['Chain A', 'Chain B', 'Chain C', 'Chain D']

for i, f_name in enumerate(files['rmsf']):
    data = read_xvg(f_name)
    if data is not None:
        res = data[:, 0]
        split_indices = np.where(np.diff(res) < 0)[0] + 1
        chains = np.split(data, split_indices)

        for c_idx, chain_data in enumerate(chains):
            if c_idx < 4:
                ax_chains[c_idx].plot(
                    chain_data[:, 0],
                    chain_data[:, 1],
                    color=colors[i],
                    label=labels[i],
                    lw=1.2,
                    alpha=0.85
                )

# Style individual chain RMSF panels
for c_idx, ax in enumerate(ax_chains):
    ax.set_title(f'E{c_idx+1}. RMSF - {chain_names[c_idx]}', fontsize=11, fontweight='bold', loc='left')
    ax.set_xlabel('Residue Index', fontsize=9)
    ax.set_ylabel('RMSF (nm)', fontsize=9)
    ax.legend(loc='upper right', frameon=True, fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.5)

plt.suptitle('Complex System Dynamics (TAN Ligand)', fontsize=16, fontweight='bold', y=0.995)

# Save figure
output_image = 'tan_complex_dynamics_separated.png'
plt.savefig(output_image, bbox_inches='tight')
print(f"Success! System dynamics saved to '{output_image}'.")