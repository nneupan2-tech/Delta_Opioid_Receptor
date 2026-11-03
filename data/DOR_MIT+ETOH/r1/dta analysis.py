import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

# Styling setup
plt.rcParams.update({'font.size': 10, 'figure.dpi': 300, 'font.sans-serif': 'Arial'})

def read_xvg(filename):
    """Parses GROMACS .xvg files while skipping header comments."""
    if not os.path.exists(filename):
        return None
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                parts = line.strip().split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
    return pd.DataFrame({'X': x, 'Y': y})

def plot_discontinuous(ax, df, label, color):
    """Plots data without drawing diagonal lines over index gaps."""
    if df is None or df.empty:
        return
    
    diffs = df['X'].diff()
    split_indices = df.index[diffs > 1].tolist()
    
    blocks = []
    prev_idx = 0
    for idx in split_indices:
        blocks.append(df.loc[prev_idx:idx-1])
        prev_idx = idx
    blocks.append(df.loc[prev_idx:])
    
    first_block = True
    for block in blocks:
        if not block.empty:
            ax.plot(
                block['X'], 
                block['Y'], 
                label=label if first_block else "", 
                color=color, 
                alpha=0.85, 
                linewidth=1.2
            )
            first_block = False

replicates = ['r1', 'r2', 'r3']
colors = {'r1': '#1f77b4', 'r2': '#ff7f0e', 'r3': '#2ca02c'}

# =========================================================
# 1. 2x2 SUMMARY GRID (MITETOH)
# =========================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle('GROMACS Simulation Analysis Summary (Mitetoh)', fontsize=15, fontweight='bold')

# Panel 1: Protein RMSD
ax_rmsd = axes[0, 0]
for rep in replicates:
    df = read_xvg(f'rmsd_mitetoh_{rep}.xvg')
    if df is not None:
        ax_rmsd.plot(df['X'], df['Y'], label=f'Replicate {rep.upper()}', color=colors[rep], alpha=0.85)
ax_rmsd.set_xlabel('Time (ns)', fontweight='bold')
ax_rmsd.set_ylabel('RMSD (nm)', fontweight='bold')
ax_rmsd.set_title('Protein RMSD', fontweight='bold')
ax_rmsd.legend(loc='best')
ax_rmsd.grid(True, linestyle='--', alpha=0.5)

# Panel 2: RMSF (CHAIN 1 ONLY)
ax_rmsf = axes[0, 1]
for rep in replicates:
    df = read_xvg(f'rmsf_mitetoh_{rep}.xvg')
    if df is not None:
        chain1_length = len(df) // 4
        chain1_df = df.iloc[0:chain1_length].reset_index(drop=True)
        plot_discontinuous(ax_rmsf, chain1_df, label=f'Replicate {rep.upper()}', color=colors[rep])
ax_rmsf.set_xlabel('Residue Index', fontweight='bold')
ax_rmsf.set_ylabel('RMSF (nm)', fontweight='bold')
ax_rmsf.set_title('Residue RMSF (Chain 1)', fontweight='bold')
ax_rmsf.legend(loc='best')
ax_rmsf.grid(True, linestyle='--', alpha=0.5)

# Panel 3: Radius of Gyration
ax_gyrate = axes[1, 0]
for rep in replicates:
    df = read_xvg(f'gyrate_mitetoh_{rep}.xvg')
    if df is not None:
        ax_gyrate.plot(df['X'], df['Y'], label=f'Replicate {rep.upper()}', color=colors[rep], alpha=0.85)
ax_gyrate.set_xlabel('Time (ns)', fontweight='bold')
ax_gyrate.set_ylabel('Rg (nm)', fontweight='bold')
ax_gyrate.set_title('Radius of Gyration', fontweight='bold')
ax_gyrate.legend(loc='best')
ax_gyrate.grid(True, linestyle='--', alpha=0.5)

# Panel 4: Protein-Unk H-Bonds
ax_hbond = axes[1, 1]
has_hbond = False
for rep in replicates:
    df = read_xvg(f'hbond_prot_unk_mitetoh_{rep}.xvg')
    if df is not None:
        ax_hbond.plot(df['X'], df['Y'], label=f'Replicate {rep.upper()}', color=colors[rep], alpha=0.85)
        has_hbond = True

ax_hbond.set_xlabel('Time (ns)', fontweight='bold')
ax_hbond.set_ylabel('Count', fontweight='bold')
if has_hbond:
    ax_hbond.set_title('Protein-Unk H-Bonds', fontweight='bold')
    ax_hbond.legend(loc='best')
else:
    ax_hbond.set_title('Protein-Unk H-Bonds (Data Partial)', fontweight='bold')
ax_hbond.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('gromacs_analysis_summary_mitetoh.png', dpi=300)
plt.close()
print("--> Saved summary grid: gromacs_analysis_summary_mitetoh.png")

# =========================================================
# 2. STANDALONE HIGH-RES RMSF CHAIN 1 PLOT
# =========================================================
fig, ax = plt.subplots(figsize=(8, 4.5))
for rep in replicates:
    df = read_xvg(f'rmsf_mitetoh_{rep}.xvg')
    if df is not None:
        chain1_length = len(df) // 4
        chain1_df = df.iloc[0:chain1_length].reset_index(drop=True)
        plot_discontinuous(ax, chain1_df, label=f'Replicate {rep.upper()}', color=colors[rep])

ax.set_xlabel('Residue Index', fontweight='bold')
ax.set_ylabel('RMSF (nm)', fontweight='bold')
ax.set_title('Protein Chain 1 RMSF (Mitetoh)', fontweight='bold')
ax.legend(frameon=True)
ax.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('rmsf_chain1_mitetoh.png', dpi=300)
plt.close()
print("--> Saved standalone plot: rmsf_chain1_mitetoh.png")