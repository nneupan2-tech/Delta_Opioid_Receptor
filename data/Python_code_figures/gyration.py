# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# SETTINGS
# ============================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Publication-quality font with bolder axis lines
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams["xtick.major.width"] = 1.3
plt.rcParams["ytick.major.width"] = 1.3


# ============================================================
# COLORS
# ============================================================

COLOR_MAP = {
    "apo": "#1f77b4",
    "drug": "#ff7f0e",
    "ethanol": "#2ca02c"
}


# ============================================================
# LABELS
# ============================================================

LABEL_MAP = {
    "apo": "Apo",
    "apoetoh": "Apo + EtOH",
    "mit": "Mitragynine",
    "mitetoh": "Mitragynine + EtOH",
    "7oh": "7-OH Mitragynine",
    "7ohetoh": "7-OH Mitragynine + EtOH",
    "tan": "TAN-67",
    "tanetoh": "TAN-67 + EtOH"
}


# ============================================================
# SYSTEM TYPE FOR COLORS
# ============================================================

SYSTEM_TYPE = {
    "apo": "apo",
    "apoetoh": "ethanol",

    "mit": "drug",
    "mitetoh": "ethanol",

    "7oh": "drug",
    "7ohetoh": "ethanol",

    "tan": "drug",
    "tanetoh": "ethanol"
}


# ============================================================
# LOAD XVG FILE
# ============================================================

def load_gyrate(filepath):
    """
    Load GROMACS gyrate.xvg file and convert Rg from nm to Angstroms.

    Expected format:
        time    Rg_x    Rg_y    Rg_z    Rg_total

    We use the final column = total radius of gyration.
    """

    time = []
    rg = []

    with open(filepath, "r", encoding="utf-8") as f:

        for line in f:

            # Skip GROMACS comments and metadata
            if line.startswith(("@", "#")):
                continue

            parts = line.split()

            if len(parts) < 2:
                continue

            try:
                t = float(parts[0])

                # Multiply by 10.0 to convert nm to Angstroms
                if len(parts) >= 5:
                    value = float(parts[4]) * 10.0
                else:
                    value = float(parts[1]) * 10.0

            except ValueError:
                continue

            time.append(t)
            rg.append(value)

    if len(time) == 0:
        raise ValueError(
            "No numerical data found in file: {}".format(filepath)
        )

    return np.array(time), np.array(rg)


# ============================================================
# LOAD REPLICATES
# ============================================================

def load_replicates(system):

    times = []
    values = []

    for replica in [1, 2, 3]:

        filename = "gyrate_{}_r{}.xvg".format(system, replica)
        filepath = os.path.join(BASE_DIR, filename)

        if not os.path.exists(filepath):

            print(
                "WARNING: {} not found - skipping".format(filename)
            )

            continue

        t, rg = load_gyrate(filepath)

        times.append(t)
        values.append(rg)

        print(
            "Loaded: {:<30} {:>8} points".format(
                filename, len(rg)
            )
        )

    if len(values) == 0:

        raise FileNotFoundError(
            "\nNo gyrate files found for system: {}\n".format(system)
        )

    print(
        "Using {} replica(s) for {}".format(
            len(values), system
        )
    )

    return times, values


# ============================================================
# CALCULATE MEAN AND SD
# ============================================================

def calculate_statistics(times, values):

    min_length = min(len(v) for v in values)

    trimmed_values = np.array([
        v[:min_length] for v in values
    ])

    trimmed_times = times[0][:min_length]

    mean = np.mean(trimmed_values, axis=0)
    sd = np.std(trimmed_values, axis=0)

    return trimmed_times, mean, sd


# ============================================================
# PLOT
# ============================================================

def generate_gyrate_plot(
    systems,
    plot_number,
    output_name
):

    plt.figure(figsize=(8.5, 5.5))

    for system in systems:

        print("\nProcessing:", system)

        times, values = load_replicates(system)

        time, mean, sd = calculate_statistics(
            times,
            values
        )

        if time.max() > 10000:
            time_ns = time / 1000.0
        else:
            time_ns = time

        color = COLOR_MAP[
            SYSTEM_TYPE[system]
        ]

        label = LABEL_MAP[system]

        # Mean line
        plt.plot(
            time_ns,
            mean,
            linewidth=2.0,
            color=color,
            label=label
        )

        # +/- 1 SD
        plt.fill_between(
            time_ns,
            mean - sd,
            mean + sd,
            color=color,
            alpha=0.18
        )

    # ========================================================
    # AXES
    # ========================================================

    plt.xlabel(
        "Time (ns)",
        fontsize=15,
        fontweight="bold",
        labelpad=8
    )

    plt.ylabel(
        "Radius of gyration (\u00C5)",
        fontsize=15,
        fontweight="bold",
        labelpad=8
    )

    plt.xlim(
        0,
        1000
    )

    # Tick line formatting
    plt.tick_params(
        axis="both",
        which="major",
        direction="out",
        length=6,
        width=1.3
    )

    # Bold and size tick labels directly
    plt.xticks(fontsize=13, fontweight="bold")
    plt.yticks(fontsize=13, fontweight="bold")

    plt.legend(
        fontsize=11,
        frameon=False,
        loc="best"
    )

    plt.grid(
        False
    )

    plt.tight_layout()

    # Save high-resolution figure
    output_path = os.path.join(
        BASE_DIR,
        output_name
    )

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        "\nSaved:",
        output_path
    )


# ============================================================
# COMPARISONS
# ============================================================

COMPARISONS = [

    (
        ["apo", "apoetoh"],
        "1",
        "gyrate_1_apo_vs_apoetoh.png"
    ),

    (
        ["apo", "mit", "mitetoh"],
        "2",
        "gyrate_2_apo_mit_mitetoh.png"
    ),

    (
        ["apo", "7oh", "7ohetoh"],
        "3",
        "gyrate_3_apo_7oh_7ohetoh.png"
    ),

    (
        ["apo", "tan", "tanetoh"],
        "4",
        "gyrate_4_apo_tan_tanetoh.png"
    )

]


# ============================================================
# MAIN
# ============================================================

print("\n")
print("=" * 65)
print("        RADIUS OF GYRATION ANALYSIS")
print("=" * 65)
print("\n")

for systems, number, output_name in COMPARISONS:

    print("\n")
    print("=" * 65)
    print("Plot", number)
    print("=" * 65)

    generate_gyrate_plot(
        systems,
        number,
        output_name
    )


print("\n")
print("=" * 65)
print("ALL RADIUS OF GYRATION PLOTS COMPLETED")
print("=" * 65)
print("\n")