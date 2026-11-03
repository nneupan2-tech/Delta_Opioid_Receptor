import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Global Plot Styling
plt.rcParams.update({'font.size': 10, 'figure.dpi': 300, 'font.sans-serif': 'Arial'})

def read_xvg(filename):
    """Parses GROMACS .xvg files safely, ignoring comments."""
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

def plot_rmsf_discontinuous(ax, df, label, color):
    """
    Plots RMSF without drawing diagonal lines across index gaps.
    Splits data whenever the gap between consecutive residue indices > 1.
    """
    if df is None or df.empty:
        return
    
    # Calculate residue index differences
    diffs = df['X'].diff()
    # Identify indices where a jump occurs
    split_indices = df.index[diffs > 1].tolist()
    
    # Split dataframe into contiguous blocks
    blocks = []
    prev_idx = 0
    for idx in split_indices:
        blocks.append(df.loc[prev_idx:idx-1])
        prev_idx = idx
    blocks.append(df.loc[prev_idx:])
    
    # Plot each continuous segment separately (only add legend to first segment)
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
# 1. GENERATE FIXED 2x2 SUMMARY PLOT
# =========================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle('GROMACS Simulation Analysis Summary', fontsize=15, fontweight='bold')

# Standard metrics
summary_configs = [
    (axes[0, 0], 'rmsd_apoetoh_{}.xvg', 'Time (ns)', 'RMSD (nm)', 'Protein RMSD'),
    (axes[1, 0], 'gyrate_apoetoh_{}.xvg', 'Time (ns)', 'Rg (nm)', 'Radius of Gyration'),
    (axes[1, 1], 'hbond_prot_etoh_apoetoh_{}.xvg', 'Time (ns)', 'Count', 'Protein-Ethanol H-Bonds')
]

for ax, pattern, xlabel, ylabel, title in summary_configs:
    for rep in replicates:
        df = read_xvg(pattern.format(rep))
        if df is not None:
            ax.plot(df['X'], df['Y'], label=f'Replicate {rep.upper()}', color=colors[rep], alpha=0.85, linewidth=1.2)
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, linestyle='--', alpha=0.5)

# Handle RMSF panel specially using gap-splitting
ax_rmsf = axes[0, 1]
for rep in replicates:
    df = read_xvg(f'rmsf_apoetoh_{rep}.xvg')
    plot_rmsf_discontinuous(ax_rmsf, df, label=f'Replicate {rep.upper()}', color=colors[rep])

ax_rmsf.set_xlabel('Residue Index', fontweight='bold')
ax_rmsf.set_ylabel('RMSF (nm)', fontweight='bold')
ax_rmsf.set_title('Residue RMSF (All Chains)', fontweight='bold')
ax_rmsf.legend(loc='best')
ax_rmsf.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('gromacs_analysis_summary_apoetoh.png', dpi=300)
plt.close()
print("--> Saved fixed summary plot without diagonal lines!")

# =========================================================
# 2. GENERATE FIXED 4 INDIVIDUAL CHAIN RMSF PLOTS
# =========================================================
sample_df = read_xvg('rmsf_apoetoh_r1.xvg')

if sample_df is not None:
    total_residues = len(sample_df)
    chain_len = total_residues // 4
    chain_names = ['Chain A', 'Chain B', 'Chain C', 'Chain D']

    for i, name in enumerate(chain_names):
        start_idx = i * chain_len
        end_idx = (i + 1) * chain_len if i < 3 else total_residues
        
        fig, ax = plt.subplots(figsize=(8, 4.5))
        
        for rep in replicates:
            df = read_xvg(f'rmsf_apoetoh_{rep}.xvg')
            if df is not None:
                # Get chain data and reset local index for slice
                chain_df = df.iloc[start_idx:end_idx].reset_index(drop=True)
                plot_rmsf_discontinuous(ax, chain_df, label=f'Replicate {rep.upper()}', color=colors[rep])
        
        ax.set_xlabel('Residue Index', fontweight='bold')
        ax.set_ylabel('RMSF (nm)', fontweight='bold')
        ax.set_title(f'Protein {name} RMSF', fontweight='bold')
        ax.legend(frameon=True)
        ax.grid(True, linestyle='--', alpha=0.5)
        
        output_name = f'rmsf_chain_{chr(65+i)}.png'
        plt.tight_layout()
        plt.savefig(output_name, dpi=300)
        plt.close()
        print(f"--> Saved clean chain plot: {output_name}")