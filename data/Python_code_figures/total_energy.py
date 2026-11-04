# -*- coding: utf-8 -*-
"""
Per-Residue MM/PBSA Energy Decomposition Plotter
Publication-Ready Formatting for JCIM / ACS Standards (Panels b, c, e, f, h, i)
"""

import os
import glob
import math
import re
import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ==============================================================================
# 1. BALLESTEROS-WEINSTEIN (BW) MAPPING & AMINO ACID DICTIONARY
# ==============================================================================
BW_MAPPING = {
    95: "2.50",   108: "ECL1", 125: "3.29",  127: "3.31",  128: "3.32",
    129: "3.33",  131: "3.35", 132: "3.36",  135: "3.39",  192: "ECL2",
    213: "5.38",  214: "5.39", 217: "5.42",  218: "5.43",  274: "6.48",
    277: "6.51",  281: "6.55", 282: "6.56",  284: "6.58",  300: "7.35",
    303: "7.38",  304: "7.39", 308: "7.43",  310: "7.45",
}

AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C",
    "ASH": "D", "GLH": "E", "LYN": "K", "ARN": "R",
}

INCLUDE_RESNUM = True
ENERGY_THRESHOLD = 0.5  # Filter threshold (|kcal/mol|)
COLOR_BAR = "#2A6F97"     # Deep Steel Blue

# ==============================================================================
# 2. PUBLICATION GRAPHICS CONFIGURATION (JCIM / ACS STANDARDS - FONT HANDLING)
# ==============================================================================
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

mpl.rcParams["font.sans-serif"] = [selected_font, 'Liberation Sans', 'Arial', 'DejaVu Sans', 'sans-serif']
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["font.size"] = 9.0
mpl.rcParams["axes.labelsize"] = 9.0
mpl.rcParams["axes.titlesize"] = 9.5
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.titleweight"] = "bold"
mpl.rcParams["xtick.labelsize"] = 8.0
mpl.rcParams["ytick.labelsize"] = 8.0
mpl.rcParams["axes.edgecolor"] = "#222222"
mpl.rcParams["axes.linewidth"] = 0.9
mpl.rcParams["pdf.fonttype"] = 42

# Mathtext configuration
mpl.rcParams["mathtext.fontset"] = "dejavusans"

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================
def clean_residue_name(res_str):
    """Formats 'A:ASP:128' -> 'D128$^{3.32}$'."""
    s = str(res_str).strip()
    s = re.sub(r"^(?:chain\s*)?[A-Za-z0-9][:_\-\s]+", "", s, flags=re.IGNORECASE).strip()

    code, num = None, None
    match = re.search(r"([A-Za-z]{3})[:_\-\s]*(\d+)", s)
    if match:
        code, num = match.groups()
        code = code.upper()
    else:
        match_alt = re.search(r"([A-Za-z]+)[:_\-\s]*(\d+)", s)
        if match_alt:
            c, n = match_alt.groups()
            code = c.upper()
            num = n

    if code and num:
        res_num = int(num)
        aa_single = AA_3TO1.get(code, code[0] if len(code) == 1 else code)

        if res_num in BW_MAPPING:
            bw = BW_MAPPING[res_num]
            res_prefix = f"{aa_single}{num}" if INCLUDE_RESNUM else f"{aa_single}"
            return rf"{res_prefix}$^{{{bw}}}$"

        return f"{aa_single}{num}"

    return s


def is_protein_residue(res_str):
    res_str_clean = str(res_str).strip().upper()
    excluded = ["UNK", "TAN", "TAN67", "MIT", "7OH", "ETOH", "LIG", "WAT", "HOH", "SOL"]
    if any(ex in res_str_clean for ex in excluded):
        return False
    codes = re.findall(r"[A-Z]{3}", res_str_clean)
    return any(code in AA_3TO1 for code in codes)


def parse_mmpbsa_deltas_tdc(filepath):
    rows = []
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        in_deltas, in_tdc, header_found = False, False, False

        for line in lines:
            line_str = line.strip()
            if "DELTAS:" in line_str:
                in_deltas = True
                continue

            if in_deltas:
                if "Total Decomposition Contribution (TDC)" in line_str:
                    in_tdc = True
                    continue

                if in_tdc:
                    if line_str.startswith("Frame #") or ("Residue" in line_str and "Internal" in line_str):
                        header_found = True
                        continue

                    if header_found:
                        if line_str.startswith(("---", "===", "Sidechain", "Backbone")) or not line_str:
                            if len(rows) > 0:
                                break
                            continue

                        parts = [p.strip() for p in line_str.split(",")]
                        if len(parts) >= 8:
                            try:
                                res = parts[1]
                                if is_protein_residue(res):
                                    rows.append({
                                        "Residue": clean_residue_name(res),
                                        "total": float(parts[7].split("+/-")[0].strip()),
                                    })
                                continue
                            except (ValueError, IndexError):
                                pass

                        pipe_parts = [p.strip() for p in line_str.split("|")]
                        if len(pipe_parts) >= 7:
                            res = pipe_parts[0]
                            if any(k in res.lower() for k in ["residue", "total", "deltas"]) or not is_protein_residue(res):
                                continue
                            try:
                                rows.append({
                                    "Residue": clean_residue_name(res),
                                    "total": float(pipe_parts[6].split("+/-")[0].strip()),
                                })
                                continue
                            except (ValueError, IndexError):
                                pass
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return pd.DataFrame(rows)


def extract_res_num(res_str):
    nums = re.findall(r"\d+", str(res_str))
    return int(nums[0]) if nums else 0


def load_and_aggregate_system(pattern):
    files = sorted(glob.glob(pattern)) or sorted(glob.glob(pattern.lower()))
    if not files:
        return None, f"No files match: {pattern}"

    df_list = [df for f in files if not (df := parse_mmpbsa_deltas_tdc(f)).empty]
    if not df_list:
        return None, "Parsing returned empty dataset"

    combined = pd.concat(df_list, axis=0)

    summary = (
        combined.groupby("Residue")
        .agg(
            total_mean=("total", "mean"),
            total_sd=("total", "std"),
            n_obs=("total", "count"),
        )
        .reset_index()
    )

    summary["total_sem"] = summary["total_sd"] / np.sqrt(summary["n_obs"].clip(lower=1))
    summary["res_num"] = summary["Residue"].apply(extract_res_num)
    summary = summary.sort_values("res_num").drop(columns=["res_num"])

    filtered = summary[summary["total_mean"].abs() >= ENERGY_THRESHOLD].reset_index(drop=True)
    return filtered, None


# ==============================================================================
# 4. SYSTEM DEFINITIONS & PLOTTING
# ==============================================================================
SYSTEMS = [
    {"panel": "b", "title": "MIT",         "pattern": "FINAL_DECOMP_MMPBSA_mit_r*.dat"},
    {"panel": "c", "title": "MIT (+EtOH)",  "pattern": "FINAL_DECOMP_MMPBSA_mitetoh_r*.dat"},
    {"panel": "e", "title": "7OH",         "pattern": "FINAL_DECOMP_MMPBSA_7oh_r*.dat"},
    {"panel": "f", "title": "7OH (+EtOH)",  "pattern": "FINAL_DECOMP_MMPBSA_7ohetoh_r*.dat"},
    {"panel": "h", "title": "TAN",         "pattern": "FINAL_DECOMP_MMPBSA_tan_r*.dat"},
    {"panel": "i", "title": "TAN (+EtOH)",  "pattern": "FINAL_DECOMP_MMPBSA_tanetoh_r*.dat"},
]

def generate_figure():
    n_systems = len(SYSTEMS)
    ncols = 2
    nrows = math.ceil(n_systems / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(8.2, 3.4 * nrows), dpi=600)
    axes = np.array(axes).flatten()

    bar_height = 0.58
    clean_error_style = dict(capsize=1.5, elinewidth=0.8, ecolor="#333333", alpha=0.90)

    for idx, sys in enumerate(SYSTEMS):
        ax = axes[idx]
        data, err_msg = load_and_aggregate_system(sys["pattern"])

        # Panel Identifier: Font Size 17, Bold
        ax.text(
            -0.16, 1.05, sys["panel"],
            transform=ax.transAxes,
            fontsize=17, fontweight="bold",
            va="bottom", ha="right"
        )

        if data is None or data.empty:
            ax.set_title(f"{sys['title']}\n(No Protein Residue Data)", fontsize=8.5, fontweight="bold", color="#666666")
            ax.set_xlim(-8.5, 5.0)
            ax.axvline(0, color="#222222", linestyle="-", linewidth=0.8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            continue

        y = np.arange(len(data["Residue"]))

        ax.barh(
            y, data["total_mean"],
            height=bar_height,
            xerr=data["total_sem"].fillna(0),
            color=COLOR_BAR,
            edgecolor="#111111",
            linewidth=0.5,
            zorder=3,
            error_kw=clean_error_style,
        )

        ax.set_yticks(y)
        ax.set_yticklabels(data["Residue"], fontsize=8.0, fontweight="bold")
        ax.invert_yaxis()

        # Synchronized X-axis across all subplots
        ax.set_xlim(-8.5, 5.0)
        ax.axvline(0, color="#222222", linestyle="-", linewidth=0.8, zorder=4)
        ax.grid(True, linestyle=":", alpha=0.45, axis="x", color="#A0A0A0")

        ax.set_title(sys["title"], pad=6, fontweight="bold", fontsize=9.5)
        ax.tick_params(axis="both", which="major", labelsize=8.0, labelbottom=True)

        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight("bold")

        if idx >= (nrows - 1) * ncols:
            ax.set_xlabel(r"Total $\Delta E$ (kcal/mol)", fontweight="bold", fontsize=9.0)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for j in range(n_systems, len(axes)):
        fig.delaxes(axes[j])

    plt.subplots_adjust(wspace=0.32, hspace=0.38, left=0.14, right=0.96, top=0.95, bottom=0.07)

    plt.savefig("Figure_DELTAS_Publication_Refined.png", dpi=600, bbox_inches="tight")
    plt.savefig("Figure_DELTAS_Publication_Refined.pdf", bbox_inches="tight")
    plt.close()
    print("[SUCCESS] Publication-ready figure (Panels b, c, e, f, h, i) generated at 600 DPI.")

if __name__ == "__main__":
    generate_figure()