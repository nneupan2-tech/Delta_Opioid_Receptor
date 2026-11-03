import numpy as np
import matplotlib.pyplot as plt

def read_xvg(filename):
    """Parses GROMACS .xvg files ignoring header comments."""
    data = []
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                cols = line.strip().split()
                if len(cols) >= 2:
                    data.append([float(c) for c in cols[:2]])
    return np.array(data)

# Set up figure grid (2x2)
fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=300)
plt.subplots_adjust(hspace=0.35, wspace=0.25)

# 1. RMSD Plot
try:
    rmsd = read_xvg('rmsd.xvg')
    axes[0, 0].plot(rmsd[:, 0], rmsd[:, 1], color='#1f77b4', lw=1)
    axes[0, 0].set_title('Protein RMSD', fontsize=12, fontweight='bold')
    axes[0, 0].set_xlabel('Time (ns)')
    axes[0, 0].set_ylabel('RMSD (nm)')
    axes[0, 0].grid(True, linestyle='--', alpha=0.5)
except Exception as e:
    axes[0, 0].set_title('RMSD (File not found)')

# 2. RMSF Plot (Splits multi-chain trajectories)
try:
    rmsf = read_xvg('rmsf.xvg')
    res = rmsf[:, 0]
    
    # Detect chain resets where residue index drops back to start
    split_indices = np.where(np.diff(res) < 0)[0] + 1
    chains = np.split(rmsf, split_indices)
    
    for idx, chain in enumerate(chains):
        axes[0, 1].plot(chain[:, 0], chain[:, 1], lw=1.2, label=f'Chain {chr(65+idx)}')
        
    axes[0, 1].set_title('Protein RMSF (Separated Chains)', fontsize=12, fontweight='bold')
    axes[0, 1].set_xlabel('Residue Number')
    axes[0, 1].set_ylabel('RMSF (nm)')
    axes[0, 1].legend(loc='upper right', fontsize=8)
    axes[0, 1].grid(True, linestyle='--', alpha=0.5)
except Exception as e:
    axes[0, 1].set_title('RMSF (File not found)')

# 3. Radius of Gyration (Rg)
try:
    rg = read_xvg('gyrate.xvg')
    axes[1, 0].plot(rg[:, 0], rg[:, 1], color='#2ca02c', lw=1)
    axes[1, 0].set_title('Radius of Gyration (Rg)', fontsize=12, fontweight='bold')
    axes[1, 0].set_xlabel('Time (ns)')
    axes[1, 0].set_ylabel('Rg (nm)')
    axes[1, 0].grid(True, linestyle='--', alpha=0.5)
except Exception as e:
    axes[1, 0].set_title('Rg (File not found)')

# 4. Protein-Ethanol (ETOH) H-Bonds
try:
    hb_etoh = read_xvg('hbond_prot_etoh.xvg')
    axes[1, 1].plot(hb_etoh[:, 0], hb_etoh[:, 1], color='#9467bd', lw=1)
    axes[1, 1].set_title('Protein - Ethanol (ETOH) H-Bonds', fontsize=12, fontweight='bold')
    axes[1, 1].set_xlabel('Time (ns)')
    axes[1, 1].set_ylabel('Number of H-Bonds')
    axes[1, 1].grid(True, linestyle='--', alpha=0.5)
except Exception as e:
    axes[1, 1].set_title('Ethanol H-Bonds (File not found)')

plt.suptitle('GROMACS Simulation Analysis Summary', fontsize=16, fontweight='bold', y=0.98)
plt.savefig('gromacs_analysis_summary.png', bbox_inches='tight')
print("Successfully saved analysis plot to 'gromacs_analysis_summary.png'!")