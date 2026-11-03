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

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle('GROMACS Simulation Analysis Summary (Tanetoh)', fontsize=15, fontweight='bold')

# 1. Protein RMSD
ax_rmsd = axes[0, 0]
for rep in replicates:
    df = read_xvg(f'rmsd_tanetoh_{rep}.xvg')
    if df is not None:
        ax_rmsd.plot(df['X'], df['Y'], label=f'Replicate {rep.upper()}', color=colors[rep], alpha=0.85)
ax_rmsd.set_xlabel('Time (ns)', fontweight='bold')
ax_rmsd.set_ylabel('RMSD (nm)', fontweight='bold')
ax_rmsd.set_title('Protein RMSD', fontweight='bold')
ax_rmsd.legend(loc='best')
ax_rmsd.grid(True, linestyle='--', alpha=0.5)

# 2. Residue RMSF (CHAIN 1 ONLY)
ax_rmsf = axes[0, 1]
for rep in replicates:
    df = read_xvg(f'rmsf_tanetoh_{rep}.xvg')
    if df is not None:
        # Slice only the first quarter of residues (Chain 1)
        total_residues = len(df)
        chain1_length = total_residues // 4
        chain1_df = df.iloc[0:chain1_length].reset_index(drop=True)
        
        plot_discontinuous(ax_rmsf, chain1_df, label=f'Replicate {rep.upper()}', color=colors[rep])

ax_rmsf.set_xlabel('Residue Index', fontweight='bold')
ax_rmsf.set_ylabel('RMSF (nm)', fontweight='bold')
ax_rmsf.set_title('Residue RMSF (Chain 1 Only)', fontweight='bold')
ax_rmsf.legend(loc='best')
ax_rmsf.grid(True, linestyle='--', alpha=0.5)

# 3. Radius of Gyration
ax_gyrate = axes[1, 0]
for rep in replicates:
    df = read_xvg(f'gyrate_tanetoh_{rep}.xvg')
    if df is not None:
        ax_gyrate.plot(df['X'], df['Y'], label=f'Replicate {rep.upper()}', color=colors[rep], alpha=0.85)
ax_gyrate.set_xlabel('Time (ns)', fontweight='bold')
ax_gyrate.set_ylabel('Rg (nm)', fontweight='bold')
ax_gyrate.set_title('Radius of Gyration', fontweight='bold')
ax_gyrate.legend(loc='best')
ax_gyrate.grid(True, linestyle='--', alpha=0.5)

# 4. Empty placeholder or optional 4th metric
ax_empty = axes[1, 1]
ax_empty.axis('off')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('gromacs_analysis_summary_tanetoh.png', dpi=300)
plt.show()
print("Successfully generated plot: gromacs_analysis_summary_tanetoh.png")