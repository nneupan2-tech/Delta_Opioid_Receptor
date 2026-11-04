# -*- coding: utf-8 -*-
"""
Publication-Quality Combined Figure Generator for JCIM (ACS Standard)
Panel a: Backbone RMSD Time Series (0-1000 ns) [2x2]
Panel b: Final 250 ns Equilibrated RMSD Statistics
Panel c: dOR Backbone RMSF Profiles (Chain 1) [2x2]
"""

import os
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy import stats

# ============================================================
# MATPLOTLIB GLOBAL STYLING (JCIM / ACS STANDARDS)
# ============================================================
# Suppress Matplotlib's font manager log warnings
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# Check available fonts and pick best available sans-serif
system_fonts = set(f.name for f in fm.fontManager.ttflist)
font_family = 'Liberation Sans' if 'Liberation Sans' in system_fonts else ('Arial' if 'Arial' in system_fonts else 'DejaVu Sans')

mpl.rcParams['font.sans-serif'] = [font_family, 'Liberation Sans', 'Arial', 'DejaVu Sans', 'sans-serif']
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

# Color Palette
COLOR_MAP = {
    'apo': '#1f77b4',       # Blue (Apo Baseline)
    'mit': '#ff7f0e',       # Orange (Drug Only)
    '7oh': '#ff7f0e',       # Orange (Drug Only)
    'tan': '#ff7f0e',       # Orange (Drug Only)
    'apoetoh': '#2ca02c',   # Green (Ethanol Conditions)
    'mitetoh': '#2ca02c',   # Green (Ethanol Conditions)
    '7ohetoh': '#2ca02c',   # Green (Ethanol Conditions)
    'tanetoh': '#2ca02c',   # Green (Ethanol Conditions)
}

LABEL_MAP = {
    'apo': 'Apo',
    'apoetoh': 'Apo + EtOH',
    'mit': 'MIT',
    'mitetoh': 'MIT + EtOH',
    '7oh': '7OH',
    '7ohetoh': '7OH + EtOH',
    'tan': 'TAN',
    'tanetoh': 'TAN + EtOH'
}

# File Mappings
rmsd_files = {
    'apo': ['dOR_backbone_RMSD_apo_r1.xvg', 'dOR_backbone_RMSD_apo_r2.xvg', 'dOR_backbone_RMSD_apo_r3.xvg'],
    'apoetoh': ['rmsd_dOR_backbone_apoetoh_r1.xvg' if os.path.exists('rmsd_dOR_backbone_apoetoh_r1.xvg') else 'dOR_backbone_RMSD_apoetoh_r1.xvg', 'dOR_backbone_RMSD_apoetoh_r2.xvg', 'dOR_backbone_RMSD_apoetoh_r3.xvg'],
    'mit': ['dOR_backbone_RMSD_mit_r1.xvg', 'dOR_backbone_RMSD_mit_r2.xvg', 'dOR_backbone_RMSD_mit_r3.xvg'],
    'mitetoh': ['dOR_backbone_RMSD_mitetoh_r1.xvg', 'dOR_backbone_RMSD_mitetoh_r2.xvg', 'dOR_backbone_RMSD_mitetoh_r3.xvg'],
    '7oh': ['dOR_backbone_RMSD_7oh_r1.xvg', 'dOR_backbone_RMSD_7oh_r2.xvg', 'dOR_backbone_RMSD_7oh_r3.xvg'],
    '7ohetoh': ['dOR_backbone_RMSD_7ohetoh_r1.xvg', 'dOR_backbone_RMSD_7ohetoh_r2.xvg', 'dOR_backbone_RMSD_7ohetoh_r3.xvg'],
    'tan': ['dOR_backbone_RMSD_tan_r1.xvg', 'dOR_backbone_RMSD_tan_r2.xvg', 'dOR_backbone_RMSD_tan_r3.xvg'],
    'tanetoh': ['dOR_backbone_RMSD_tanetoh_r1.xvg', 'dOR_backbone_RMSD_tanetoh_r2.xvg', 'dOR_backbone_RMSD_tanetoh_r3.xvg'],
}

rmsf_files = {
    'apo': ['rmsf_apo_r1.xvg', 'rmsf_apo_r2.xvg', 'rmsf_apo_r3.xvg'],
    'apoetoh': ['rmsf_apoetoh_r1.xvg', 'rmsf_apoetoh_r2.xvg', 'rmsf_apoetoh_r3.xvg'],
    'mit': ['rmsf_mit_r1.xvg', 'rmsf_mit_r2.xvg', 'rmsf_mit_r3.xvg'],
    'mitetoh': ['rmsf_mitetoh_r1.xvg', 'rmsf_mitetoh_r2.xvg', 'rmsf_mitetoh_r3.xvg'],
    '7oh': ['rmsf_7oh_r1.xvg', 'rmsf_7oh_r2.xvg', 'rmsf_7oh_r3.xvg'],
    '7ohetoh': ['rmsf_7ohetoh_r1.xvg', 'rmsf_7ohetoh_r2.xvg', 'rmsf_7ohetoh_r3.xvg'],
    'tan': ['rmsf_tan_r1.xvg', 'rmsf_tan_r2.xvg', 'rmsf_tan_r3.xvg'],
    'tanetoh': ['rmsf_tanetoh_r1.xvg', 'rmsf_tanetoh_r2.xvg', 'rmsf_tanetoh_r3.xvg']
}

# ============================================================
# DATA LOADING FUNCTIONS
# ============================================================
def load_xvg_rmsd(filepath):
    if not os.path.exists(filepath):
        return np.array([]), np.array([])
    time, rmsd = [], []
    unit_is_ps = False
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('@'):
                if 'xaxis' in line and '(ps)' in line:
                    unit_is_ps = True
                continue
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                rmsd.append(float(parts[1]) * 10.0)  # Convert nm to Angstroms
    time, rmsd = np.array(time), np.array(rmsd)
    if unit_is_ps or (len(time) > 0 and np.max(time) > 2000):
        time = time / 1000.0
    return time, rmsd

def load_all_chains_rmsf(filepath):
    if not os.path.exists(filepath):
        return []
    all_blocks, current_residues, current_rmsf = [], [], []
    previous_residue = None
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('@') or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                residue = float(parts[0])
                value = float(parts[1]) * 10.0  # Convert nm to Angstroms
            except ValueError:
                continue
            if previous_residue is not None and residue < previous_residue and len(current_residues) > 0:
                all_blocks.append((np.array(current_residues), np.array(current_rmsf)))
                current_residues, current_rmsf = [], []
            current_residues.append(residue)
            current_rmsf.append(value)
            previous_residue = residue
    if len(current_residues) > 0:
        all_blocks.append((np.array(current_residues), np.array(current_rmsf)))
    return all_blocks

def load_dor_rmsf(filepath):
    blocks = load_all_chains_rmsf(filepath)
    if not blocks:
        return np.array([]), np.array([])
    return blocks[0]  # Return Chain 1 (dOR)

def collect_equilibrated_rmsd(system_name, window_ns=250.0):
    files = rmsd_files.get(system_name, [])
    rep_means = []
    for f in files:
        t_ns, r = load_xvg_rmsd(f)
        if len(t_ns) > 0:
            max_time = np.max(t_ns)
            window_rmsd = r[t_ns >= (max_time - window_ns)]
            if len(window_rmsd) > 0:
                rep_means.append(np.mean(window_rmsd))
    return np.array(rep_means)

def draw_significance_bracket(ax, x1, x2, y, h, p_val, text_offset=0.03):
    sig_symbol = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color='black', lw=1.0)
    ax.text((x1 + x2) * 0.5, y + h + text_offset, sig_symbol,
            ha='center', va='bottom', fontsize=8.5, fontweight='bold')

# ============================================================
# MASTER PLOTTING FUNCTION
# ============================================================
def generate_combined_figure(output_filename='Figure2_Combined_dOR_Dynamics.png'):
    fig = plt.figure(figsize=(7.2, 11.2), dpi=600)
    
    gs_main = gridspec.GridSpec(3, 1, height_ratios=[2.2, 1.25, 2.2], figure=fig, 
                                hspace=0.38, left=0.10, right=0.98, top=0.95, bottom=0.05)
    
    # --- PANEL A: RMSD TIME SERIES ---
    gs_a = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_main[0], wspace=0.15, hspace=0.28)
    comparisons_a = [
        (['apo', 'apoetoh'], 'Apo vs Apo + EtOH'),
        (['apo', 'mit', 'mitetoh'], 'Apo vs MIT vs MIT + EtOH'),
        (['apo', '7oh', '7ohetoh'], 'Apo vs 7OH vs 7OH + EtOH'),
        (['apo', 'tan', 'tanetoh'], 'Apo vs TAN vs TAN + EtOH')
    ]
    
    for idx, (systems, title) in enumerate(comparisons_a):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs_a[row, col])
        
        if idx == 0:
            ax.text(-0.16, 1.10, 'a', transform=ax.transAxes, fontsize=15, fontweight='bold', va='bottom', ha='right')
            
        for sys in systems:
            files = rmsd_files.get(sys, [])
            valid_files = [f for f in files if os.path.exists(f)]
            if not valid_files:
                continue
            time_series, rmsd_series = [], []
            for f in valid_files:
                t, r = load_xvg_rmsd(f)
                if len(t) > 0:
                    time_series.append(t)
                    rmsd_series.append(r)
            if rmsd_series:
                min_len = min(len(r) for r in rmsd_series)
                t_axis = time_series[0][:min_len]
                stacked = np.vstack([r[:min_len] for r in rmsd_series])
                mean_rmsd = np.mean(stacked, axis=0)
                std_rmsd = np.std(stacked, axis=0)
                c = COLOR_MAP.get(sys, '#333333')
                
                ax.plot(t_axis, mean_rmsd, label=LABEL_MAP.get(sys, sys), lw=1.2, color=c)
                ax.fill_between(t_axis, mean_rmsd - std_rmsd, mean_rmsd + std_rmsd, alpha=0.18, color=c)
        
        ax.set_ylim(0, 4.8)
        ax.set_xlim(0, 1000)
        ax.set_title(title, fontsize=9.0, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.25)
        ax.tick_params(labelsize=8.5)
        
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight('bold')
        
        ax.legend(loc='upper left', fontsize=7.0, frameon=False, handlelength=1.5, handletextpad=0.4)
        
        if col == 0:
            ax.set_ylabel(r'Backbone RMSD ($\mathrm{\AA}$)', fontsize=9.0, fontweight='bold')
        else:
            ax.set_yticklabels([])
            
        if row == 1:
            ax.set_xlabel('Time (ns)', fontsize=9.0, fontweight='bold')
        else:
            ax.set_xticklabels([])

    # --- PANEL B: EQUILIBRATED STATISTICAL BAR CHART ---
    ax_b = fig.add_subplot(gs_main[1])
    ax_b.text(-0.08, 1.10, 'b', transform=ax_b.transAxes, fontsize=15, fontweight='bold', va='bottom', ha='right')
    
    systems_order = ['apo', 'apoetoh', 'mit', 'mitetoh', '7oh', '7ohetoh', 'tan', 'tanetoh']
    replicate_data = {sys: collect_equilibrated_rmsd(sys) for sys in systems_order}
    group_means = [np.mean(replicate_data[s]) if len(replicate_data[s]) > 0 else np.nan for s in systems_order]
    group_stds = [np.std(replicate_data[s], ddof=1) if len(replicate_data[s]) > 1 else 0.0 for s in systems_order]
    
    x_pos = np.arange(len(systems_order))
    colors = [COLOR_MAP.get(s, '#333333') for s in systems_order]
    
    ax_b.bar(x_pos, group_means, yerr=group_stds, align='center', alpha=0.75,
             ecolor='black', capsize=3, color=colors, edgecolor='black', lw=1.0, width=0.45, zorder=2)
    
    np.random.seed(42)
    for i, sys in enumerate(systems_order):
        reps = replicate_data[sys]
        if len(reps) > 0:
            jitter = np.random.normal(0, 0.02, size=len(reps))
            ax_b.scatter(np.full_like(reps, i) + jitter, reps, color='black', edgecolor='white', lw=0.8, s=35, zorder=6)

    all_reps = [v for v in replicate_data.values() if len(v) > 0]
    y_max = np.max(np.concatenate(all_reps)) if all_reps else 3.8
    
    etoh_comps = [('apo', 'apoetoh'), ('mit', 'mitetoh'), ('7oh', '7ohetoh'), ('tan', 'tanetoh')]
    y_lvl1 = y_max + 0.50
    bracket_h = 0.12
    
    for s1, s2 in etoh_comps:
        idx1, idx2 = systems_order.index(s1), systems_order.index(s2)
        r1, r2 = replicate_data[s1], replicate_data[s2]
        if len(r1) > 1 and len(r2) > 1:
            _, p_val = stats.ttest_ind(r1, r2, equal_var=False)
            draw_significance_bracket(ax_b, idx1, idx2, y_lvl1, bracket_h, p_val)

    apo_drug_comps = [('apo', 'mit'), ('apo', '7oh'), ('apo', 'tan')]
    y_curr = y_lvl1 + 0.65
    for s1, s2 in apo_drug_comps:
        idx1, idx2 = systems_order.index(s1), systems_order.index(s2)
        r1, r2 = replicate_data[s1], replicate_data[s2]
        if len(r1) > 1 and len(r2) > 1:
            _, p_val = stats.ttest_ind(r1, r2, equal_var=False)
            draw_significance_bracket(ax_b, idx1, idx2, y_curr, bracket_h, p_val)
            y_curr += 0.65

    x_labels = ['APO', 'APO+ETOH', 'MIT', 'MIT+ETOH', '7OH', '7OH+ETOH', 'TAN', 'TAN+ETOH']
    ax_b.set_ylabel(r'Final 250 ns RMSD ($\mathrm{\AA}$)', fontsize=9.0, fontweight='bold')
    ax_b.set_title('Equilibrated Dynamics (750-1000 ns): Replicate Means & Grand Means', fontsize=9.0, fontweight='bold', pad=8)
    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(x_labels, fontsize=8.5, fontweight='bold')
    ax_b.set_ylim(0, y_curr + 0.4)
    ax_b.tick_params(axis='both', labelsize=8.5)
    for tick in ax_b.get_xticklabels() + ax_b.get_yticklabels():
        tick.set_fontweight('bold')
    ax_b.grid(True, axis='y', linestyle='--', alpha=0.3, zorder=1)
    
    legend_patches = [
        Patch(facecolor='#1f77b4', edgecolor='black', alpha=0.75, label='Apo Baseline'),
        Patch(facecolor='#ff7f0e', edgecolor='black', alpha=0.75, label='Drug Conditions'),
        Patch(facecolor='#2ca02c', edgecolor='black', alpha=0.75, label='Ethanol Conditions')
    ]
    ax_b.legend(handles=legend_patches, loc='upper right', fontsize=7.5, frameon=False, handlelength=1.2)

    # --- PANEL C: dOR BACKBONE RMSF PROFILES ---
    gs_c = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_main[2], wspace=0.15, hspace=0.28)
    
    for idx, (systems, title) in enumerate(comparisons_a):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs_c[row, col])
        
        if idx == 0:
            ax.text(-0.16, 1.10, 'c', transform=ax.transAxes, fontsize=15, fontweight='bold', va='bottom', ha='right')
            
        for sys in systems:
            files = rmsf_files.get(sys, [])
            valid_files = [f for f in files if os.path.exists(f)]
            if not valid_files:
                continue
            res_series, rmsf_series = [], []
            for f in valid_files:
                res, rmsf = load_dor_rmsf(f)
                if len(res) > 0:
                    res_series.append(res)
                    rmsf_series.append(rmsf)
            if rmsf_series:
                min_len = min(len(r) for r in rmsf_series)
                res_axis = res_series[0][:min_len]
                stacked = np.vstack([r[:min_len] for r in rmsf_series])
                mean_rmsf = np.mean(stacked, axis=0)
                std_rmsf = np.std(stacked, axis=0)
                c = COLOR_MAP.get(sys, '#333333')
                
                ax.plot(res_axis, mean_rmsf, label=LABEL_MAP.get(sys, sys), lw=1.2, color=c)
                ax.fill_between(res_axis, mean_rmsf - std_rmsf, mean_rmsf + std_rmsf, alpha=0.18, color=c)
        
        ax.set_ylim(0, 9.0)
        ax.set_title(title, fontsize=9.0, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.25)
        ax.tick_params(labelsize=8.5)
        
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight('bold')
            
        ax.legend(loc='upper left', fontsize=7.0, frameon=False, handlelength=1.5, handletextpad=0.4)
        
        if col == 0:
            ax.set_ylabel(r'$\mathbf{\delta}$OR RMSF ($\mathrm{\AA}$)', fontsize=9.0, fontweight='bold')
        else:
            ax.set_yticklabels([])
            
        if row == 1:
            ax.set_xlabel('Residue Number', fontsize=9.0, fontweight='bold')
        else:
            ax.set_xticklabels([])

    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Figure saved cleanly to '{output_filename}' at 600 DPI.")

if __name__ == '__main__':
    generate_combined_figure('Figure2_Combined_dOR_Dynamics.png')