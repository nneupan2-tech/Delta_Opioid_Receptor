# -*- coding: utf-8 -*-
"""
Publication-Quality Composite Figure Generator for JCIM / ACS Standards
Panels a-f formatted with Arial 9 pt and enhanced 15 pt Panel Labels.
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
from scipy import stats
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ============================================================
# 1. GLOBAL STYLING (ACS STANDARDS - FONT HANDLING)
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
mpl.rcParams['font.size'] = 9.0
mpl.rcParams['axes.titlesize'] = 9.0
mpl.rcParams['axes.labelsize'] = 9.0
mpl.rcParams['xtick.labelsize'] = 9.0
mpl.rcParams['ytick.labelsize'] = 9.0
mpl.rcParams['legend.fontsize'] = 8.0

mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['xtick.major.size'] = 3.5
mpl.rcParams['ytick.major.size'] = 3.5

# Mathtext configuration
mpl.rcParams['mathtext.fontset'] = 'dejavusans'

# Standardized ACS Color Palette
COLOR_APO     = "#1f77b4"  # Blue
COLOR_DRUG    = "#ff7f0e"  # Orange
COLOR_ALCOHOL = "#2ca02c"  # Green

COLOR_MAP = {
    'apo': COLOR_APO,     'apoetoh': COLOR_ALCOHOL,
    'mit': COLOR_DRUG,    'mitetoh': COLOR_ALCOHOL,
    '7oh': COLOR_DRUG,    '7ohetoh': COLOR_ALCOHOL,
    'tan': COLOR_DRUG,    'tanetoh': COLOR_ALCOHOL
}

FILE_GROUPS = {
    'apo':     [f'D253_displacement_over_time_apo_r{i}.csv' for i in [1, 2, 3]],
    'apoetoh': [f'D253_displacement_over_time_apoetoh_r{i}.csv' for i in [1, 2, 3]],
    'mit':     [f'D253_displacement_over_time_mit_r{i}.csv' for i in [1, 2, 3]],
    'mitetoh': [f'D253_displacement_over_time_mitetoh_r{i}.csv' for i in [1, 2, 3]],
    '7oh':     [f'D253_displacement_over_time_7oh_r{i}.csv' for i in [1, 2, 3]],
    '7ohetoh': [f'D253_displacement_over_time_7ohetoh_r{i}.csv' for i in [1, 2, 3]],
    'tan':     [f'D253_displacement_over_time_tan_r{i}.csv' for i in [1, 2, 3]],
    'tanetoh': [f'D253_displacement_over_time_tanetoh_r{i}.csv' for i in [1, 2, 3]],
}

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================
def load_and_filter_data(filepath, t_start=750, t_end=1000):
    if not os.path.exists(filepath):
        return np.array([])
    df = pd.read_csv(filepath)
    filtered = df[(df['time_ns'] >= t_start) & (df['time_ns'] <= t_end)]
    return filtered['displacement_A'].values

def collect_group_data(sys_name, t_start=750, t_end=1000):
    files = FILE_GROUPS.get(sys_name, [])
    rep_means, pooled = [], []
    for f in files:
        disp = load_and_filter_data(f, t_start, t_end)
        if len(disp) > 0:
            rep_means.append(np.mean(disp))
            pooled.extend(disp)
    return np.array(pooled), np.array(rep_means)

def add_stat_bracket(ax, x1, x2, p_val, y_pos, bar_height=0.25):
    if p_val < 0.001:
        text = '***'
    elif p_val < 0.01:
        text = '**'
    elif p_val < 0.05:
        text = '*'
    else:
        text = 'ns'
        
    ax.plot([x1, x1, x2, x2], [y_pos, y_pos + bar_height, y_pos + bar_height, y_pos], lw=1.0, c='black')
    ax.text((x1 + x2) * 0.5, y_pos + bar_height + 0.10, text, ha='center', va='bottom', color='black', fontsize=9.0, fontweight='bold')

# ============================================================
# 3. MASTER PLOTTING FUNCTION
# ============================================================
def generate_master_figure(output_filename='master_d253_composite_arial9.png', t_start=750, t_end=1000):
    fig = plt.figure(figsize=(7.2, 9.8), dpi=600)
    gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1.30], hspace=0.38, wspace=0.28)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    # Panel Identifiers (Bold 15 pt for prominent visual structure)
    panel_axes = [('a', ax_a), ('b', ax_b), ('c', ax_c), ('d', ax_d), ('e', ax_e), ('f', ax_f)]
    for label, ax in panel_axes:
        ax.text(
            -0.12, 1.08, label,
            transform=ax.transAxes,
            fontsize=15.0,
            fontweight='bold',
            va='bottom',
            ha='right'
        )

    # --------------------------------------------------------
    # PANELS A-D: TIME-SERIES TRAJECTORIES
    # --------------------------------------------------------
    ts_configs = [
        (ax_a, "Apo vs Apo+EtOH", [('apo', 'APO', COLOR_APO), ('apoetoh', 'APO+EtOH', COLOR_ALCOHOL)]),
        (ax_b, "Apo vs MIT vs MIT+EtOH", [('apo', 'APO', COLOR_APO), ('mit', 'MIT', COLOR_DRUG), ('mitetoh', 'MIT+EtOH', COLOR_ALCOHOL)]),
        (ax_c, "Apo vs 7OH vs 7OH+EtOH", [('apo', 'APO', COLOR_APO), ('7oh', '7OH', COLOR_DRUG), ('7ohetoh', '7OH+EtOH', COLOR_ALCOHOL)]),
        (ax_d, "Apo vs TAN vs TAN+EtOH", [('apo', 'APO', COLOR_APO), ('tan', 'TAN', COLOR_DRUG), ('tanetoh', 'TAN+EtOH', COLOR_ALCOHOL)])
    ]

    for ax, title, series in ts_configs:
        for sys_key, label, color in series:
            files = FILE_GROUPS[sys_key]
            all_r = []
            ref_t = None
            for f in files:
                if os.path.exists(f):
                    df = pd.read_csv(f)
                    all_r.append(df['displacement_A'].values)
                    if ref_t is None:
                        ref_t = df['time_ns'].values
            if all_r:
                min_l = min(len(r) for r in all_r)
                stacked = np.vstack([r[:min_l] for r in all_r])
                mean_v = pd.Series(np.mean(stacked, axis=0)).rolling(50, min_periods=1, center=True).mean()
                std_v = pd.Series(np.std(stacked, axis=0)).rolling(50, min_periods=1, center=True).mean()
                
                ax.plot(ref_t[:min_l], mean_v, label=label, color=color, lw=1.2)
                ax.fill_between(ref_t[:min_l], mean_v - std_v, mean_v + std_v, color=color, alpha=0.15)

        ax.set_title(title, fontsize=9.0, fontweight='bold', pad=4)
        ax.set_xlim(0, 1000)
        ax.set_ylim(7.5, 16.5)
        ax.grid(True, linestyle='--', alpha=0.25)
        ax.tick_params(axis='both', labelsize=9.0)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight('bold')
            
        ax.legend(loc='upper right', frameon=False, fontsize=8.0)

    # Axis Labels for Panels a-d (Bold 9 pt)
    ax_a.set_ylabel(r'D253 Displacement ($\mathrm{\AA}$)', fontsize=9.0, fontweight='bold')
    ax_c.set_ylabel(r'D253 Displacement ($\mathrm{\AA}$)', fontsize=9.0, fontweight='bold')
    ax_c.set_xlabel('Time (ns)', fontsize=9.0, fontweight='bold')
    ax_d.set_xlabel('Time (ns)', fontsize=9.0, fontweight='bold')

    # --------------------------------------------------------
    # PANEL E: GPCR SCHEMATIC PLACEHOLDER
    # --------------------------------------------------------
    ax_e.set_title("GPCR Conformational Transition", fontsize=9.0, fontweight='bold', pad=4)
    ax_e.text(0.5, 0.5, "Insert GPCR Structural Model\n(INACTIVE vs ACTIVE State)", 
              ha='center', va='center', fontsize=9.0, fontweight='bold', color='gray')
    ax_e.set_xticks([])
    ax_e.set_yticks([])
    for spine in ax_e.spines.values():
        spine.set_linestyle('--')
        spine.set_color('gray')

    # --------------------------------------------------------
    # PANEL F: EQUILIBRATED BOXPLOT (750-1000 ns)
    # --------------------------------------------------------
    all_systems = ['apo', 'apoetoh', 'mit', 'mitetoh', '7oh', '7ohetoh', 'tan', 'tanetoh']
    system_labels = ['APO', 'APO\n+EtOH', 'MIT', 'MIT\n+EtOH', '7OH', '7OH\n+EtOH', 'TAN', 'TAN\n+EtOH']

    data_to_plot, rep_means_list, sys_idx_map = [], [], {}
    for idx, sys in enumerate(all_systems, start=1):
        pooled, r_means = collect_group_data(sys, t_start, t_end)
        data_to_plot.append(pooled if len(pooled) > 0 else np.array([12.0]))
        rep_means_list.append(r_means)
        sys_idx_map[sys] = idx

    bp = ax_f.boxplot(
        data_to_plot, patch_artist=True, labels=system_labels, widths=0.36, showmeans=True,
        meanprops={"marker": "D", "markeredgecolor": "black", "markerfacecolor": "white", "markersize": 4.0},
        medianprops={"color": "black", "linewidth": 1.2},
        flierprops={"marker": "o", "markersize": 1.5, "alpha": 0.15}
    )

    colors = [COLOR_MAP[s] for s in all_systems]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.70)
        patch.set_edgecolor('black')

    for idx, r_means in enumerate(rep_means_list, start=1):
        if len(r_means) > 0:
            x_jit = np.random.normal(idx, 0.03, size=len(r_means))
            ax_f.scatter(x_jit, r_means, color='black', edgecolor='white', linewidth=0.6, zorder=4, s=22)

    # Pairwise statistical tests & brackets
    pairwise_tests = [('apo', 'apoetoh'), ('mit', 'mitetoh'), ('7oh', '7ohetoh'), ('tan', 'tanetoh')]
    all_max_y = []
    for sys_a, sys_b in pairwise_tests:
        x1, x2 = sys_idx_map[sys_a], sys_idx_map[sys_b]
        _, rep_a = collect_group_data(sys_a, t_start, t_end)
        _, rep_b = collect_group_data(sys_b, t_start, t_end)
        
        if len(rep_a) > 1 and len(rep_b) > 1:
            _, p_val = stats.ttest_ind(rep_a, rep_b, equal_var=False)
            pair_max = max(np.max(data_to_plot[x1-1]), np.max(data_to_plot[x2-1]))
            all_max_y.append(pair_max)
            add_stat_bracket(ax_f, x1, x2, p_val, y_pos=pair_max + 0.40, bar_height=0.20)

    # EXPAND HEADROOM TO SEPARATE BRACKETS & LEGEND COMPLETELY
    highest_y = max(all_max_y) if all_max_y else 15.5
    ax_f.set_ylim(7.5, highest_y + 4.8)

    # Vertical Condition Dividers
    ax_f.axvline(x=2.5, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)
    ax_f.axvline(x=4.5, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)
    ax_f.axvline(x=6.5, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)

    ax_f.set_ylabel(r'D253 Displacement ($\mathrm{\AA}$)', fontsize=9.0, fontweight='bold')
    ax_f.set_title(f"Equilibrated State ({t_start}-{t_end} ns)", fontsize=9.0, fontweight='bold', pad=4)
    ax_f.tick_params(axis='both', which='major', labelsize=9.0)
    for tick in ax_f.get_xticklabels() + ax_f.get_yticklabels():
        tick.set_fontweight('bold')
    ax_f.grid(True, axis='y', linestyle='--', alpha=0.25)

    # Borderless, Two-Column Legend Positioned Safely at Upper Right
    legend_elements = [
        Patch(facecolor=COLOR_APO, edgecolor='black', alpha=0.7, label='Apo'),
        Patch(facecolor=COLOR_DRUG, edgecolor='black', alpha=0.7, label='Drug'),
        Patch(facecolor=COLOR_ALCOHOL, edgecolor='black', alpha=0.7, label='+EtOH'),
        Line2D([0], [0], marker='D', color='w', markeredgecolor='black', markerfacecolor='white', markersize=4.5, label='Grand Mean'),
        Line2D([0], [0], marker='o', color='w', markeredgecolor='white', markerfacecolor='black', markersize=4.5, label='Replicate Mean')
    ]
    ax_f.legend(handles=legend_elements, loc='upper right', fontsize=8.0, frameon=False, ncol=2, columnspacing=0.8)

    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Master composite figure saved to '{output_filename}' at 600 DPI.")

if __name__ == '__main__':
    generate_master_figure()