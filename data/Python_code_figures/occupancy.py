# -*- coding: utf-8 -*-
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Presentation style configuration
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'figure.dpi': 300
})

# System mappings matching your filename keywords
systems = {
    '7ohetoh': '7-OH-Mitragynine',
    'mitetoh': 'Mitragynine',
    'tanetoh': 'TAN-67'
}

colors = {
    '7-OH-Mitragynine': '#D95F02',  # Orange/Red
    'Mitragynine': '#1B9E77',       # Teal
    'TAN-67': '#7570B3'             # Purple
}

def load_metric(prefix, is_occupancy=False):
    """Finds all replicates (r1, r2, r3) for each system and averages them."""
    data_dict = {}
    for sys_key, sys_name in systems.items():
        pattern = f"{prefix}_*{sys_key}_r*.xvg"
        files = sorted(glob.glob(pattern))
        
        if not files:
            print(f"Warning: No files matched pattern '{pattern}'")
            continue
            
        reps = []
        time = None
        for f in files:
            raw = np.loadtxt(f, comments=['#', '@'])
            time = raw[:, 0] / 1000.0  # Convert ps to ns
            vals = raw[:, 1]
            if is_occupancy:
                vals = vals / 9.0  # Convert 9 atoms to ~1 ethanol molecule
            reps.append(vals)
            
        rep_matrix = np.array(reps)
        data_dict[sys_name] = {
            'time': time,
            'mean': np.mean(rep_matrix, axis=0),
            'std': np.std(rep_matrix, axis=0)
        }
    return data_dict

def smooth_series(data_array, window=50):
    """Applies a moving average to smooth thermal noise (50 frames = 5 ns)."""
    return pd.Series(data_array).rolling(window=window, min_periods=1).mean().to_numpy()

# 1. LOAD DATA
print("Loading metrics across all replicates...")
mindist = load_metric("mindist_ethanol_drug")
numcont = load_metric("numcont_ethanol_drug")
occupancy = load_metric("ethanol_occupancy", is_occupancy=True)

# 2. PLOT 1: MINIMUM DISTANCE
plt.figure(figsize=(9, 4.5))
for sys_name, data in mindist.items():
    smoothed_mean = smooth_series(data['mean'])
    plt.plot(data['time'], smoothed_mean, label=sys_name, color=colors[sys_name], lw=2)
    plt.fill_between(data['time'], smooth_series(data['mean'] - data['std']), 
                     smooth_series(data['mean'] + data['std']), color=colors[sys_name], alpha=0.18)

plt.axhline(0.4, color='gray', linestyle='--', alpha=0.7, label='Cutoff (0.4 nm)')
plt.ylabel('Min Distance (nm)')
plt.xlabel('Time (ns)')
plt.title('Minimum Distance Between Drug and Ethanol (5 ns Moving Avg)', fontweight='bold')
plt.legend(loc='upper right', frameon=True)
plt.tight_layout()
plt.savefig('minimum_distance.png', dpi=300)
print("Saved figure: minimum_distance.png")
plt.close()

# 3. PLOT 2: CONTACT COUNT
plt.figure(figsize=(9, 4.5))
for sys_name, data in numcont.items():
    smoothed_mean = smooth_series(data['mean'])
    plt.plot(data['time'], smoothed_mean, label=sys_name, color=colors[sys_name], lw=2)
    plt.fill_between(data['time'], smooth_series(data['mean'] - data['std']), 
                     smooth_series(data['mean'] + data['std']), color=colors[sys_name], alpha=0.18)

plt.ylabel('Contact Count')
plt.xlabel('Time (ns)')
plt.title('Atomic Contacts Within Solvation Shell (<0.4 nm)', fontweight='bold')
plt.legend(loc='upper right', frameon=True)
plt.tight_layout()
plt.savefig('contact_count.png', dpi=300)
print("Saved figure: contact_count.png")
plt.close()

# 4. PLOT 3: ETHANOL OCCUPANCY
plt.figure(figsize=(9, 4.5))
for sys_name, data in occupancy.items():
    smoothed_mean = smooth_series(data['mean'])
    plt.plot(data['time'], smoothed_mean, label=sys_name, color=colors[sys_name], lw=2)
    plt.fill_between(data['time'], smooth_series(data['mean'] - data['std']), 
                     smooth_series(data['mean'] + data['std']), color=colors[sys_name], alpha=0.18)

plt.ylabel('Ethanol Molecules')
plt.xlabel('Time (ns)')
plt.ylim(bottom=0)
plt.title('Solvation Shell Ethanol Occupancy', fontweight='bold')
plt.legend(loc='upper right', frameon=True)
plt.tight_layout()
plt.savefig('ethanol_occupancy.png', dpi=300)
print("Saved figure: ethanol_occupancy.png")
plt.close()