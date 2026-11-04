# -*- coding: utf-8 -*-
"""
Publication-Quality Figure Generator for JCIM / ACS Standards
Panel b: Gai alpha5-Helix Relative RMSD Dynamics (0-1000 ns) [2x2 Grid]
"""

import os
import glob
import re
import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ============================================================
# 1. MATPLOTLIB GLOBAL STYLING (JCIM / ACS STANDARDS - FONT HANDLING)
# ============================================================
# Suppress Matplotlib's font manager log warnings
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# Automatically register valid TrueType fonts in user's ~/.local/share/fonts directory
user_fonts_dir = os.path.expanduser('~/.local/share/fonts')
if os.path.exists(user_fonts_dir):
    for font_file in os.listdir(user_fonts_dir):
        font_path = os.path.join(user_fonts_dir, font_file)
        if font_file.lower().endswith(('.ttf', '.otf')) and os.path.getsize(font_path) > 0:
            try:
                fm.fontManager.addfont(font_path)
            except Exception:
                pass

# Select Liberation Sans (Arial metrical equivalent) if available, fallback to DejaVu Sans
available_fonts = set(f.name for f in fm.fontManager.ttflist)
selected_font = 'Liberation Sans' if 'Liberation Sans' in available_fonts else ('Arial' if 'Arial' in available_fonts else 'DejaVu Sans')

mpl.rcParams['font.sans-serif'] = [selected_font, 'Liberation Sans', 'Arial', 'DejaVu Sans', 'sans-serif']
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['xtick.major.size'] = 3.5
mpl.rcParams['ytick.major.size'] = 3.5

# Mathtext configuration
mpl.rcParams['mathtext.fontset'] = 'dejavusans'

DATA_DIR = "."
SMOOTH_WINDOW = 50  # Centered moving average window (frames)

# Color Palette
COLOR_APO     = "#1f77b4"  # Blue
COLOR_DRUG    = "#ff7f0e"  # Orange
COLOR_ALCOHOL = "#2ca02c"  # Green

# Clean Subplot Configurations (NO "System:" and NO main figure title)
SUBPLOT_CONFIGS = [
    ("Apo vs Apo+EtOH", [("APO", "APO", COLOR_APO, "-"), ("APOETOH", "APO+EtOH", COLOR_ALCOHOL, "-")]),
    ("Apo vs 7OH vs 7OH+EtOH", [("APO", "APO", COLOR_APO, "-"), ("7OH", "7OH", COLOR_DRUG, "-"), ("7OHETOH", "7OH+EtOH", COLOR_ALCOHOL, "-")]),
    ("Apo vs MIT vs MIT+EtOH", [("APO", "APO", COLOR_APO, "-"), ("MIT", "MIT", COLOR_DRUG, "-"), ("MITETOH", "MIT+EtOH", COLOR_ALCOHOL, "-")]),
    ("Apo vs TAN vs TAN+EtOH", [("APO", "APO", COLOR_APO, "-"), ("TAN", "TAN", COLOR_DRUG, "-"), ("TANETOH", "TAN+EtOH", COLOR_ALCOHOL, "-")])
]

# ============================================================
# 2. DATA EXTRACTION & HELPER FUNCTIONS
# ============================================================
def read_xvg(file_path):
    """Parses GROMACS .xvg file, converts time (ps to ns), and distance (nm to Angstrom)."""
    data = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("@", "#")) and line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    time_ns = float(parts[0]) / 1000.0  # ps to ns
                    dist = float(parts[1])
                    data.append([time_ns, dist])
    if not data:
        return np.array([]), np.array([])
    
    df = pd.DataFrame(data, columns=["Time_ns", "Distance"])
    
    # Auto-detect nm vs Angstroms
    if df["Distance"].max() < 2.0:
        df["Distance_A"] = df["Distance"] * 10.0
    else:
        df["Distance_A"] = df["Distance"]
        
    return df["Time_ns"].values, df["Distance_A"].values

def moving_average_pandas(data, window_size):
    """Applies centered rolling mean."""
    if window_size <= 1:
        return data
    return pd.Series(data).rolling(window=window_size, min_periods=1, center=True).mean().to_numpy()

def find_condition_files(condition_key):
    """Locates all replicate files for a given condition key in a case-insensitive manner."""
    all_files = glob.glob(os.path.join(DATA_DIR, "alpha5_relative_rmsd_*.xvg"))
    matched_files = []
    
    for filepath in all_files:
        filename = os.path.basename(filepath)
        match = re.search(r"alpha5_relative_rmsd_(.*)_r\d+\.xvg", filename, re.IGNORECASE)
        if match:
            cond = match.group(1).upper().replace("-", "").replace("_", "")
            target = condition_key.upper().replace("-", "").replace("_", "")
            if cond == target:
                matched_files.append(filepath)
                
    return sorted(matched_files)

def process_condition_replicates(condition_key):
    """Finds replicates, extracts time series, computes mean & std."""
    files = find_condition_files(condition_key)
    if not files:
        print(f"[WARNING] No files found matching condition: '{condition_key}'")
        return None, None, None
        
    all_distances, ref_time = [], None
    for f in files:
        t, r = read_xvg(f)
        if len(t) > 0:
            all_distances.append(r)
            if ref_time is None:
                ref_time = t
                
    if not all_distances:
        return None, None, None
        
    min_len = min(len(d) for d in all_distances)
    ref_time = ref_time[:min_len]
    stacked = np.vstack([d[:min_len] for d in all_distances])
    
    mean_raw = np.mean(stacked, axis=0)
    std_raw = np.std(stacked, axis=0)
    
    return ref_time, moving_average_pandas(mean_raw, SMOOTH_WINDOW), moving_average_pandas(std_raw, SMOOTH_WINDOW)

# ============================================================
# 3. MAIN PLOTTING FUNCTION
# ============================================================
def generate_clean_alpha5_grid(output_filename="alpha5_rmsd_clean_grid_600dpi.png"):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), sharex=True, sharey=True, dpi=600)
    axes_flat = axes.flatten()
    
    for idx, (title, series_list) in enumerate(SUBPLOT_CONFIGS):
        ax = axes_flat[idx]
        row, col = divmod(idx, 2)
        
        # Threshold reference lines
        ax.axhline(y=2.5, color='gray', linestyle='--', linewidth=0.9, alpha=0.75, zorder=1)
        ax.axhline(y=4.0, color='red', linestyle=':', linewidth=0.9, alpha=0.75, zorder=1)
        
        for sys_key, label, color, linestyle in series_list:
            t, mean_r, std_r = process_condition_replicates(sys_key)
            if t is not None:
                ax.plot(t, mean_r, label=label, color=color, linestyle=linestyle, lw=1.2, zorder=3)
                ax.fill_between(t, mean_r - std_r, mean_r + std_r, color=color, alpha=0.15, zorder=2)
                
        ax.set_ylim(0.0, 5.8)
        ax.set_xlim(0, 1000)
        
        # Set standardized titles (No "System:" prefixes)
        ax.set_title(title, fontsize=9.0, fontweight='bold', pad=5)
        ax.grid(True, linestyle="--", alpha=0.25)
        
        # Tick parameters enlarged and bolded
        ax.tick_params(axis='both', which='major', labelsize=8.0)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight('bold')
        
        # Clean borderless legend
        ax.legend(loc="upper right", frameon=False, fontsize=7.0, handlelength=1.4, handletextpad=0.4)
        
        # Enlarged, Bold Axis Labels
        if col == 0:
            ax.set_ylabel(r"G$\mathrm{\alpha_i}$ $\mathrm{\alpha5}$-Helix RMSD ($\mathrm{\AA}$)", fontsize=9.5, fontweight='bold')
        if row == 1:
            ax.set_xlabel("Time (ns)", fontsize=9.5, fontweight='bold')

    # Add single text callouts for threshold reference lines on top-right panel
    axes[0, 1].text(980, 2.6, r'Bound (< 2.5 $\mathrm{\AA}$)', fontsize=6.5, color='gray', ha='right', va='bottom', fontweight='bold')
    axes[0, 1].text(980, 4.1, r'Uncoupled (> 4.0 $\mathrm{\AA}$)', fontsize=6.5, color='red', ha='right', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Clean grid saved without system headers to '{output_filename}' at 600 DPI.")

if __name__ == "__main__":
    generate_clean_alpha5_grid()