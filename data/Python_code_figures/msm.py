# -*- coding: utf-8 -*-
"""
Markov State Model (MSM) Analysis, Validation & Visualization Pipeline
Targeted for High-Impact Physical Chemistry & Biophysics Publications (JCIM Format)
"""

import glob
import logging
import os
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Deeptime imports
from deeptime.clustering import KMeans
from deeptime.decomposition import TICA
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.plots import plot_ck_test, plot_implied_timescales
from deeptime.util.validation import implied_timescales
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# 0. Global Publication Aesthetics & Physical Units Configuration (JCIM Standard)
# -----------------------------------------------------------------------------
FRAME_TO_NS = 0.1  # Conversion factor: 1 frame = 0.1 ns (100 ps/frame)
MSM_LAGTIME_FRAMES = 15
MSM_LAGTIME_NS = MSM_LAGTIME_FRAMES * FRAME_TO_NS  # 1.5 ns

# Suppress Matplotlib's font manager log warnings
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# Automatically register valid TrueType fonts in user's ~/.local/share/fonts directory
user_fonts_dir = os.path.expanduser("~/.local/share/fonts")
if os.path.exists(user_fonts_dir):
    for font_file in os.listdir(user_fonts_dir):
        font_path = os.path.join(user_fonts_dir, font_file)
        if font_file.lower().endswith((".ttf", ".otf")) and os.path.getsize(font_path) > 0:
            try:
                fm.fontManager.addfont(font_path)
            except Exception:
                pass

# Select Liberation Sans (Arial metrical equivalent) if available, fallback to DejaVu Sans
available_fonts = set(f.name for f in fm.fontManager.ttflist)
selected_font = (
    "Liberation Sans"
    if "Liberation Sans" in available_fonts
    else ("Arial" if "Arial" in available_fonts else "DejaVu Sans")
)

# Global Matplotlib settings for publication figure formatting
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [selected_font, "Liberation Sans", "Arial", "DejaVu Sans", "sans-serif"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.labelweight": "bold",
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 9,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.linewidth": 1.2,
        "lines.linewidth": 1.5,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "mathtext.fontset": "dejavusans",
    }
)

# -----------------------------------------------------------------------------
# 1. Configuration & Data Loading
# -----------------------------------------------------------------------------
DATA_DIR = "/project/das/Netra/delta/Markvo_state_analysis"
CONDITIONS = [
    "apo",
    "apoetoh",
    "7oh",
    "7ohetoh",
    "mit",
    "mitetoh",
    "tan",
    "tanetoh",
]
REPLICAS = ["r1", "r2", "r3"]

DISPLAY_TITLES = {
    "apo": "APO",
    "apoetoh": "APO+EtOH",
    "7oh": "7OH",
    "7ohetoh": "7OH+EtOH",
    "mit": "MIT",
    "mitetoh": "MIT+EtOH",
    "tan": "TAN",
    "tanetoh": "TAN+EtOH",
}
PANEL_LABELS = ["a", "b", "c", "d", "e", "f", "g", "h"]


def read_xvg(filepath):
    """Parses GROMACS .xvg file, ignoring header comments."""
    values = []
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("@"):
                continue
            values.append(float(line.split()[1]))
    return np.array(values)


def read_csv(filepath):
    """Parses displacement .csv file, returning the last column as values."""
    df = pd.read_csv(filepath)
    return df.iloc[:, -1].values


print("--> Loading and assembling feature matrices...")
feature_trajectories = []
trajectory_keys = []

for cond in CONDITIONS:
    for rep in REPLICAS:
        tm3_tm6_file = os.path.join(
            DATA_DIR, f"tm3_tm6_dist_{cond}_{rep}.xvg"
        )
        alpha5_file = os.path.join(
            DATA_DIR, f"alpha5_relative_rmsd_{cond}_{rep}.xvg"
        )
        d253_file = os.path.join(
            DATA_DIR, f"D253_displacement_over_time_{cond}_{rep}.csv"
        )

        if not (
            os.path.exists(tm3_tm6_file)
            and os.path.exists(alpha5_file)
            and os.path.exists(d253_file)
        ):
            print(
                f"Warning: Missing files for '{cond}' replica '{rep}'. Skipping."
            )
            continue

        f1_tm3_tm6 = read_xvg(tm3_tm6_file)
        f2_alpha5 = read_xvg(alpha5_file)
        f3_d253 = read_csv(d253_file)

        min_len = min(len(f1_tm3_tm6), len(f2_alpha5), len(f3_d253))

        feat_matrix = np.column_stack(
            [f1_tm3_tm6[:min_len], f3_d253[:min_len], f2_alpha5[:min_len]]
        )

        feature_trajectories.append(feat_matrix)
        trajectory_keys.append((cond, rep))

print(f"--> Successfully loaded {len(feature_trajectories)} trajectories.")

# -----------------------------------------------------------------------------
# 2. Fit Global TICA Space & Microstates
# -----------------------------------------------------------------------------
print("--> Fitting Global TICA space across all trajectories...")
TICA_LAGTIME = 10
tica_model = (
    TICA(lagtime=TICA_LAGTIME, dim=2).fit(feature_trajectories).fetch_model()
)
tica_trajectories = [tica_model.transform(f) for f in feature_trajectories]

N_CLUSTERS = 100
print(f"--> Clustering TICA trajectories into {N_CLUSTERS} microstates...")
kmeans_model = (
    KMeans(n_clusters=N_CLUSTERS, max_iter=500)
    .fit(tica_trajectories)
    .fetch_model()
)
dtrajs = [kmeans_model.transform(t) for t in tica_trajectories]

dtrajs_by_condition = {cond: [] for cond in CONDITIONS}
for idx, (cond, rep) in enumerate(trajectory_keys):
    dtrajs_by_condition[cond].append(dtrajs[idx])

# -----------------------------------------------------------------------------
# 3. MSM Fitting & Macrostate Free Energies
# -----------------------------------------------------------------------------
results = {}
print("\n=== Fitting Condition-Specific MSMs ===")
for cond, cond_dtrajs in dtrajs_by_condition.items():
    if not cond_dtrajs:
        continue
    try:
        msm_estimator = MaximumLikelihoodMSM(lagtime=MSM_LAGTIME_FRAMES)
        msm_model = msm_estimator.fit(cond_dtrajs).fetch_model()

        pcca_model = msm_model.pcca(3)
        macro_memberships = pcca_model.memberships

        pi = msm_model.stationary_distribution
        macro_pi = np.dot(macro_memberships.T, pi)

        dg = -0.593 * np.log(macro_pi + 1e-12)
        dg -= np.min(dg)

        results[cond] = {
            "msm": msm_model,
            "pcca": pcca_model,
            "macro_pi": macro_pi,
            "dg": dg,
        }

        print(
            f"Condition {cond.upper()}: Active Microstates ="
            f" {msm_model.n_states}, Macro Pop = {macro_pi.round(3)}"
        )
    except Exception as e:
        print(f"Could not fit MSM for condition {cond}: {e}")

# -----------------------------------------------------------------------------
# 4. Publication-Grade Validation Plots (ITS & Chapman-Kolmogorov)
# -----------------------------------------------------------------------------
print("\n=== Generating Harmonized Publication-Grade Validation Plots ===")
its_lagtimes_frames = np.arange(1, 50, 4)

for cond, cond_dtrajs in dtrajs_by_condition.items():
    if not cond_dtrajs:
        continue

    # ITS Plot
    try:
        its_models = []
        for l in its_lagtimes_frames:
            try:
                counts_est = TransitionCountEstimator(
                    lagtime=l, count_mode="sliding"
                )
                c_model = counts_est.fit_fetch(cond_dtrajs)
                c_conn = c_model.submodel_largest()
                msm = MaximumLikelihoodMSM().fit_fetch(c_conn)
                its_models.append(msm)
            except Exception:
                pass

        if len(its_models) > 0:
            its_data = implied_timescales(its_models)

            fig, ax = plt.subplots(figsize=(5.5, 4))
            plot_implied_timescales(its_data, n_its=3, ax=ax)

            ax.axvline(
                MSM_LAGTIME_FRAMES,
                color="red",
                linestyle="--",
                lw=1.5,
                label=f"Chosen Lag ({MSM_LAGTIME_NS:.1f} ns)",
            )

            ax.set_yscale("log")
            ax.set_title(
                f"Implied Timescales - {DISPLAY_TITLES.get(cond, cond.upper())}",
                fontweight="bold",
                fontsize=10,
            )
            ax.set_xlabel("Lag time tau (frames)", fontsize=9, fontweight="bold")
            ax.set_ylabel("Implied timescale (frames)", fontsize=9, fontweight="bold")
            ax.tick_params(axis="both", which="major", labelsize=8)
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_fontweight("bold")
            ax.grid(True, which="both", ls=":", alpha=0.5)
            ax.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=9)

            plt.savefig(f"its_{cond}.png", dpi=300, bbox_inches="tight")
            plt.close()
    except Exception as e:
        print(f"   ITS plot failed for {cond}: {e}")

    # CK Test Plot
    try:
        ck_lags_frames = [MSM_LAGTIME_FRAMES * i for i in range(1, 6)]
        ck_models = []

        for lag in ck_lags_frames:
            try:
                counts_est = TransitionCountEstimator(
                    lagtime=lag, count_mode="sliding"
                )
                c_model = counts_est.fit_fetch(cond_dtrajs)
                c_conn = c_model.submodel_largest()
                msm = MaximumLikelihoodMSM().fit_fetch(c_conn)
                ck_models.append(msm)
            except Exception:
                pass

        if len(ck_models) > 1:
            n_sets = 3
            ck_res = ck_models[0].ck_test(ck_models, n_metastable_sets=n_sets)

            plot_ck_test(ck_res)
            fig = plt.gcf()
            fig.set_size_inches(7.5, 7.0)
            fig.legends.clear()

            all_axes = fig.get_axes()

            for idx, ax_sub in enumerate(all_axes):
                ax_sub.set_ylim(0.0, 1.0)
                ax_sub.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

                ticks = ax_sub.get_xticks()
                ax_sub.set_xticks(ticks)
                ax_sub.set_xticklabels([f"{t * FRAME_TO_NS:.1f}" for t in ticks], fontsize=8, fontweight="bold")
                ax_sub.set_xlabel("Lag time (ns)", fontsize=9, fontweight="bold")
                ax_sub.set_ylabel("Probability", fontsize=9, fontweight="bold")
                ax_sub.tick_params(axis="both", which="major", labelsize=8)
                for tick in ax_sub.get_yticklabels():
                    tick.set_fontweight("bold")
                ax_sub.grid(True, ls=":", alpha=0.4)

            top_right_ax = all_axes[n_sets - 1]
            lines = top_right_ax.get_lines()
            if len(lines) >= 2:
                lines[0].set_label("Estimate (MD)")
                lines[1].set_label("Predict (MSM)")
                top_right_ax.legend(
                    loc="upper right",
                    frameon=True,
                    framealpha=0.9,
                    fontsize=8,
                )

            plt.suptitle(
                f"Chapman-Kolmogorov Test - {DISPLAY_TITLES.get(cond, cond.upper())} (tau ="
                f" {MSM_LAGTIME_NS:.1f} ns)",
                fontsize=11,
                fontweight="bold",
                y=0.98,
            )

            plt.savefig(f"ck_test_{cond}.png", dpi=300, bbox_inches="tight")
            plt.close()
    except Exception as e:
        print(f"   CK test failed for {cond}: {e}")

# -----------------------------------------------------------------------------
# 5. Comparative 2D Free Energy Map Plot with Macrostate Overlay
# -----------------------------------------------------------------------------
print(
    "\n--> Generating publication 2D Free Energy Map with Macrostate overlays..."
)
fig = plt.figure(figsize=(14, 6.5), dpi=300)
gs = fig.add_gridspec(
    2, 5, width_ratios=[1, 1, 1, 1, 0.08], wspace=0.25, hspace=0.30
)

axes = []
for r in range(2):
    for c in range(4):
        axes.append(fig.add_subplot(gs[r, c]))

cax = fig.add_subplot(gs[:, 4])

macro_colors = ["#E63946", "#F4A261", "#2A9D8F"]  # Red, Orange, Teal
macro_names = ["M0 (Uncoupled)", "M1 (Pre-coupled)", "M2 (Active)"]

cf = None
for idx, cond in enumerate(CONDITIONS):
    ax = axes[idx]
    cond_indices = [i for i, (c, r) in enumerate(trajectory_keys) if c == cond]

    if not cond_indices or cond not in results:
        continue

    cond_tica = np.vstack([tica_trajectories[i] for i in cond_indices])
    counts, xedges, yedges = np.histogram2d(
        cond_tica[:, 0], cond_tica[:, 1], bins=40
    )
    prob = counts / np.sum(counts)
    free_energy = -0.593 * np.log(prob + 1e-6)
    free_energy -= np.min(free_energy)
    free_energy = np.clip(free_energy, 0, 5)

    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    cf = ax.contourf(
        X, Y, free_energy.T, levels=np.linspace(0, 5, 21), cmap="viridis_r"
    )

    msm = results[cond]["msm"]
    pcca = results[cond]["pcca"]
    macro_pi = results[cond]["macro_pi"]
    pi = msm.stationary_distribution
    active_states = msm.state_symbols()

    active_tica_centers = kmeans_model.cluster_centers[active_states]
    weighted_memberships = pi[:, None] * pcca.memberships
    macro_tica_centers = (
        np.dot(weighted_memberships.T, active_tica_centers) / macro_pi[:, None]
    )

    for m_idx in range(3):
        x_m, y_m = macro_tica_centers[m_idx]
        ax.scatter(
            x_m,
            y_m,
            c=macro_colors[m_idx],
            s=120,
            edgecolors="black",
            linewidth=1.0,
            zorder=5,
        )
        ax.annotate(
            f"M{m_idx}",
            (x_m, y_m),
            color="white",
            fontsize=7,
            fontweight="bold",
            ha="center",
            va="center",
            zorder=6,
        )

    # Panel labels 'a' - 'h': FONT SIZE 17 BOLD
    ax.text(
        -0.14,
        1.08,
        PANEL_LABELS[idx],
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        va="bottom",
        ha="right",
    )

    # Subplot Title
    ax.set_title(
        DISPLAY_TITLES.get(cond, cond.upper()),
        fontsize=10,
        fontweight="bold",
        pad=4,
    )

    # Bold tick formatting and axis labels
    ax.tick_params(axis="both", which="major", labelsize=8)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")

    if idx >= 4:
        ax.set_xlabel("TICA 1", fontsize=9, fontweight="bold")
    if idx % 4 == 0:
        ax.set_ylabel("TICA 2", fontsize=9, fontweight="bold")

if cf is not None:
    cbar = fig.colorbar(cf, cax=cax, ticks=np.arange(0, 6, 1))
    cbar.set_label(
        "Relative Free Energy (kcal/mol)", fontsize=9, fontweight="bold"
    )
    cbar.ax.tick_params(labelsize=8)
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight("bold")

legend_elements = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        label=macro_names[i],
        markerfacecolor=macro_colors[i],
        markeredgecolor="k",
        markersize=8,
    )
    for i in range(3)
]
fig.legend(
    handles=legend_elements,
    loc="upper center",
    bbox_to_anchor=(0.46, 0.98),
    ncol=3,
    fontsize=9.5,
    frameon=False,
)

plt.suptitle(
    "Delta Opioid Receptor Free Energy Landscapes & Macrostate Centers",
    fontsize=11,
    fontweight="bold",
    y=1.02,
)
plt.savefig("msm_free_energy_landscapes.png", dpi=300, bbox_inches="tight")
plt.close()

print(
    "--> Process Complete! All plots generated to journal-quality standards."
)