# -*- coding: utf-8 -*-
"""
Publication-Quality Combined Figure Generator for JCIM / ACS Standards
TM3-TM6 Intracellular Cavity Opening Dynamics
Panel b: Time-Resolved Distance Trajectories (0-1000 ns) [2x2 Grid]
Panel c: Equilibrated Cavity Opening Statistics (750-1000 ns) [Bar Plot + Welch's t-tests]
"""

import os
import glob
import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats

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
mpl.rcParams['axes.labelweight'] = 'bold'
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['xtick.major.size'] = 3.5
mpl.rcParams['ytick.major.size'] = 3.5

# Mathtext configuration
mpl.rcParams['mathtext.fontset'] = 'dejavusans'

DATA_DIR = "."
SMOOTH_WINDOW = 50      # Trajectory moving average window (frames)
T_START = 750.0         # Start of equilibrated window (ns)
T_END = 1000.0          # End of trajectory (ns)

# Color Palette
COLOR_APO     = "#1f77b4"  # Blue
COLOR_DRUG    = "#ff7f0e"  # Orange
COLOR_ALCOHOL = "#2ca02c"  # Green

SYSTEMS = ["apo", "apoetoh", "7oh", "7ohetoh", "mit", "mitetoh", "tan", "tanetoh"]

LABELS = {
    "apo": "APO", "apoetoh": "APO+EtOH",
    "7oh": "7OH", "7ohetoh": "7OH+EtOH",
    "mit": "MIT", "mitetoh": "MIT+EtOH",
    "tan": "TAN", "tanetoh": "TAN+EtOH"
}

COLORS = {
    "apo": COLOR_APO, "apoetoh": COLOR_ALCOHOL,
    "7oh": COLOR_DRUG, "7ohetoh": COLOR_ALCOHOL,
    "mit": COLOR_DRUG, "mitetoh": COLOR_ALCOHOL,
    "tan": COLOR_DRUG, "tanetoh": COLOR_ALCOHOL
}

SUBPLOT_CONFIGS = [
    ("Apo vs Apo+EtOH", [("apo", "APO", COLOR_APO, "-"), ("apoetoh", "APO+EtOH", COLOR_ALCOHOL, "-")]),
    ("Apo vs 7OH vs 7OH+EtOH", [("apo", "APO", COLOR_APO, "-"), ("7oh", "7OH", COLOR_DRUG, "-"), ("7ohetoh", "7OH+EtOH", COLOR_ALCOHOL, "-")]),
    ("Apo vs MIT vs MIT+EtOH", [("apo", "APO", COLOR_APO, "-"), ("mit", "MIT", COLOR_DRUG, "-"), ("mitetoh", "MIT+EtOH", COLOR_ALCOHOL, "-")]),
    ("Apo vs TAN vs TAN+EtOH", [("apo", "APO", COLOR_APO, "-"), ("tan", "TAN", COLOR_DRUG, "-"), ("tanetoh", "TAN+EtOH", COLOR_ALCOHOL, "-")])
]

# Strictly direct pairwise comparisons between each condition and its +EtOH counterpart
COMPARISONS = [
    ("apo", "apoetoh"),
    ("7oh", "7ohetoh"),
    ("mit", "mitetoh"),
    ("tan", "tanetoh")
]

# ============================================================
# 2. DATA EXTRACTION & HELPER FUNCTIONS
# ============================================================
def read_xvg(filepath):
    """Loads GROMACS XVG file ignoring headers."""
    data = []
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(('@', '#')):
                data.append([float(x) for x in line.split()])
    data = np.array(data)
    if len(data) == 0:
        return np.array([]), np.array([])
    return data[:, 0], data[:, 1] * 10.0  # Convert nm to Angstroms

def moving_average_pandas(data, window_size):
    """Applies centered rolling mean."""
    if window_size <= 1:
        return data
    return pd.Series(data).rolling(window=window_size, min_periods=1, center=True).mean().to_numpy()

def process_system_replicates(system_key):
    """Finds all replicates, converts units, computes mean & std."""
    pattern = os.path.join(DATA_DIR, f"tm3_tm6_dist_{system_key}_r*.xvg")
    files = sorted(glob.glob(pattern))
    if not files:
        return None, None, None
    all_distances, ref_time = [], None
    for f in files:
        time, dist = read_xvg(f)
        if len(time) > 0:
            all_distances.append(dist)
            if ref_time is None:
                ref_time = time
    if not all_distances:
        return None, None, None
    min_len = min(len(d) for d in all_distances)
    ref_time = ref_time[:min_len]
    all_distances = np.array([d[:min_len] for d in all_distances])
    
    mean_raw = np.mean(all_distances, axis=0)
    std_raw = np.std(all_distances, axis=0)
    return ref_time, moving_average_pandas(mean_raw, SMOOTH_WINDOW), moving_average_pandas(std_raw, SMOOTH_WINDOW)

def read_xvg_window(filepath, t_start, t_end):
    """Filters trajectory within time window."""
    times, dists = [], []
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(('@', '#')):
                parts = [float(x) for x in line.split()]
                if t_start <= parts[0] <= t_end:
                    times.append(parts[0])
                    dists.append(parts[1] * 10.0)  # nm to Angstroms
    return np.array(dists)

def draw_significance_bracket(ax, x1, x2, p_val, y_pos, bar_height=0.25):
    """Draws statistical bracket with p-value symbol."""
    symbol = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
    ax.plot([x1, x1, x2, x2], [y_pos, y_pos + bar_height, y_pos + bar_height, y_pos], lw=1.0, color='black')
    ax.text((x1 + x2) * 0.5, y_pos + bar_height + 0.10, symbol, ha='center', va='bottom', color='black', fontsize=9, fontweight='bold')

# ============================================================
# 3. MASTER PLOTTING FUNCTION
# ============================================================
def generate_combined_cavity_figure(output_filename='Figure_TM3_TM6_Panels_b_c.png'):
    fig = plt.figure(figsize=(7.2, 9.2), dpi=600)
    
    # Grid layout with tight panel-to-panel spacing
    gs_main = gridspec.GridSpec(2, 1, height_ratios=[2.0, 1.15], figure=fig,
                                hspace=0.24, left=0.10, right=0.98, top=0.95, bottom=0.06)
    
    # ------------------------------------------------------------
    # PANEL B: TIME-SERIES TRAJECTORIES (2x2 GRID)
    # ------------------------------------------------------------
    gs_b = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_main[0], wspace=0.14, hspace=0.20)
    
    for idx, (title, series_list) in enumerate(SUBPLOT_CONFIGS):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs_b[row, col])
        
        # Panel label 'b' at top-left of the first subplot (Font size 15 Bold)
        if idx == 0:
            ax.text(-0.16, 1.10, 'b', transform=ax.transAxes, fontsize=15, fontweight='bold', va='bottom', ha='right')
            
        for system_key, label, color, linestyle in series_list:
            t, mean_d, std_d = process_system_replicates(system_key)
            if t is not None:
                ax.plot(t, mean_d, label=label, color=color, linestyle=linestyle, lw=1.2)
                ax.fill_between(t, mean_d - std_d, mean_d + std_d, color=color, alpha=0.15)
                
        ax.set_ylim(13.5, 18.5)
        ax.set_xlim(0, 1000)
        ax.set_title(title, fontsize=9.0, fontweight='bold', pad=4)
        ax.grid(True, linestyle="--", alpha=0.25)
        ax.tick_params(axis='both', labelsize=8.5)
        
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight('bold')
            
        # Clean, borderless legend
        ax.legend(loc="upper left", frameon=False, fontsize=7.5, handlelength=1.5, handletextpad=0.4)
        
        if col == 0:
            ax.set_ylabel(r"TM3-TM6 Distance ($\mathrm{\AA}$)", fontsize=9.0, fontweight='bold')
        else:
            ax.set_yticklabels([])
            
        if row == 1:
            ax.set_xlabel("Time (ns)", fontsize=9.0, fontweight='bold')
        else:
            ax.set_xticklabels([])

    # ------------------------------------------------------------
    # PANEL C: EQUILIBRATED STATISTICAL BAR CHART
    # ------------------------------------------------------------
    ax_c = fig.add_subplot(gs_main[1])
    # Panel label 'c' at top-left (Font size 15 Bold)
    ax_c.text(-0.08, 1.10, 'c', transform=ax_c.transAxes, fontsize=15, fontweight='bold', va='bottom', ha='right')
    
    replicate_means = {}
    results = {}
    
    for sys_key in SYSTEMS:
        pattern = os.path.join(DATA_DIR, f"tm3_tm6_dist_{sys_key}_r*.xvg")
        files = sorted(glob.glob(pattern))
        rep_means = []
        for f in files:
            dists = read_xvg_window(f, T_START, T_END)
            if len(dists) > 0:
                rep_means.append(np.mean(dists))
        rep_means = np.array(rep_means)
        if len(rep_means) > 0:
            replicate_means[sys_key] = rep_means
            results[sys_key] = {
                "mean": np.mean(rep_means),
                "std": np.std(rep_means, ddof=1) if len(rep_means) > 1 else 0.0
            }

    keys_present = [k for k in SYSTEMS if k in results]
    x_map = {k: i for i, k in enumerate(keys_present)}
    x_pos = np.arange(len(keys_present))
    means = [results[k]["mean"] for k in keys_present]
    stds = [results[k]["std"] for k in keys_present]
    bar_colors = [COLORS[k] for k in keys_present]

    # Main Bars
    ax_c.bar(x_pos, means, yerr=stds, capsize=3, color=bar_colors, alpha=0.75,
             edgecolor='black', lw=1.0, width=0.45, zorder=2)

    # Replicate Scatter Overlay
    np.random.seed(42)
    for idx, k in enumerate(keys_present):
        reps = replicate_means[k]
        jitter = np.random.normal(0, 0.02, size=len(reps))
        ax_c.scatter(np.full_like(reps, idx) + jitter, reps, color='black', edgecolor='white', lw=0.8, s=30, zorder=6)

    # Section Dividers
    for v_pos in [1.5, 3.5, 5.5]:
        ax_c.axvline(x=v_pos, color='gray', linestyle=':', alpha=0.4, lw=0.8)

    # Lift brackets comfortably higher above bar tops & scatter points
    col_top = np.array([m + s for m, s in zip(means, stds)])
    y_bracket_base = np.max(col_top) + 0.80

    # Draw direct pairwise brackets for each condition vs +EtOH
    for g1, g2 in COMPARISONS:
        if g1 in replicate_means and g2 in replicate_means:
            x1, x2 = x_map[g1], x_map[g2]
            rep_a, rep_b = replicate_means[g1], replicate_means[g2]
            p_val = 0.2  # Default ns if n < 2
            if len(rep_a) > 1 and len(rep_b) > 1:
                _, p_val = stats.ttest_ind(rep_a, rep_b, equal_var=False)
            draw_significance_bracket(ax_c, x1, x2, p_val, y_pos=y_bracket_base, bar_height=0.25)

    # Axis Formatting & Headroom Expansion
    ax_c.set_xticks(x_pos)
    ax_c.set_xticklabels([LABELS[k] for k in keys_present], fontsize=8.5, fontweight='bold')
    ax_c.set_ylabel(r"TM3-TM6 Mean Distance ($\mathrm{\AA}$)", fontsize=9.0, fontweight='bold')
    ax_c.set_title(f"Equilibrated Cavity Opening ({T_START:.0f}-{T_END:.0f} ns)", fontsize=9.0, fontweight='bold', pad=8)
    ax_c.grid(True, axis='y', linestyle='--', alpha=0.3, zorder=1)
    
    # Generous y-limit (24.0 Å) to guarantee clean gap between brackets and top-right legend
    ax_c.set_ylim(0, 24.0)
    ax_c.tick_params(axis='both', labelsize=8.5)
    for tick in ax_c.get_xticklabels() + ax_c.get_yticklabels():
        tick.set_fontweight('bold')

    # Panel c Legend (2-column format placed cleanly in upper right corner above brackets)
    legend_elements = [
        Patch(facecolor=COLOR_APO, edgecolor='black', alpha=0.75, label='APO Baseline'),
        Patch(facecolor=COLOR_DRUG, edgecolor='black', alpha=0.75, label='Drug Conditions'),
        Patch(facecolor=COLOR_ALCOHOL, edgecolor='black', alpha=0.75, label='Ethanol Conditions'),
        Line2D([0], [0], marker='o', color='w', markeredgecolor='white', markerfacecolor='black', markersize=5, label='Replicate Means')
    ]
    ax_c.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.99, 0.98), fontsize=7.5, frameon=False, handlelength=1.2, ncol=2)

    # Export High-Resolution PNG (600 DPI for publication)
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Figure saved cleanly to '{output_filename}' at 600 DPI.")

if __name__ == '__main__':
    generate_combined_cavity_figure('Figure_TM3_TM6_Panels_b_c.png')