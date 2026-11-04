import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path


# ============================================================
# MM/PBSA FRAME-BY-FRAME ANALYSIS
#
# Systems:
#   Mitragynine       +/- ethanol
#   7-OH-mitragynine  +/- ethanol
#   TAN-67            +/- ethanol
#
# Replicates:
#   3 independent replicas per condition
#
# Main quantity:
#   Mean PBSA TOTAL energy across frames
#
# Ethanol effect:
#   DeltaDeltaG = Mean(+EtOH) - Mean(-EtOH)
#
# Statistics:
#   Paired t-test across matched replicas
#
# ============================================================


# ============================================================
# 1. FILE MAPPING
# ============================================================

systems = {

    "Mitragynine_-EtOH": [
        "PER_RESIDUE_ENERGIES_mit_r1.csv",
        "PER_RESIDUE_ENERGIES_mit_r2.csv",
        "PER_RESIDUE_ENERGIES_mit_r3.csv"
    ],

    "Mitragynine_+EtOH": [
        "PER_RESIDUE_ENERGIES_mitetoh_r1.csv",
        "PER_RESIDUE_ENERGIES_mitetoh_r2.csv",
        "PER_RESIDUE_ENERGIES_mitetoh_r3.csv"
    ],

    "7OH_-EtOH": [
        "PER_RESIDUE_ENERGIES_7oh_r1.csv",
        "PER_RESIDUE_ENERGIES_7oh_r2.csv",
        "PER_RESIDUE_ENERGIES_7oh_r3.csv"
    ],

    "7OH_+EtOH": [
        "PER_RESIDUE_ENERGIES_7ohetoh_r1.csv",
        "PER_RESIDUE_ENERGIES_7ohetoh_r2.csv",
        "PER_RESIDUE_ENERGIES_7ohetoh_r3.csv"
    ],

    "TAN_-EtOH": [
        "PER_RESIDUE_ENERGIES_tan_r1.csv",
        "PER_RESIDUE_ENERGIES_tan_r2.csv",
        "PER_RESIDUE_ENERGIES_tan_r3.csv"
    ],

    "TAN_+EtOH": [
        "PER_RESIDUE_ENERGIES_tanetoh_r1.csv",
        "PER_RESIDUE_ENERGIES_tanetoh_r2.csv",
        "PER_RESIDUE_ENERGIES_tanetoh_r3.csv"
    ]
}


# ============================================================
# 2. READ MMPBSA FRAME-BY-FRAME CSV
# ============================================================

def load_mmpbsa_csv(filepath):

    filepath = Path(filepath)

    if not filepath.exists():

        raise FileNotFoundError(
            "File not found: {}".format(filepath)
        )

    header_idx = None

    # --------------------------------------------------------
    # Find the actual header
    # --------------------------------------------------------

    with open(
        filepath,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for idx, line in enumerate(f):

            if (
                "Frame #" in line
                and "TOTAL" in line
            ):

                header_idx = idx
                break

    if header_idx is None:

        raise ValueError(
            "Could not find Frame#/TOTAL header in {}".format(
                filepath
            )
        )

    # --------------------------------------------------------
    # Read CSV
    # --------------------------------------------------------

    df = pd.read_csv(
        filepath,
        skiprows=header_idx
    )

    # Clean column names
    df.columns = (
        df.columns
        .astype(str)
        .str.strip()
    )

    # --------------------------------------------------------
    # Convert numeric columns
    # --------------------------------------------------------

    for col in df.columns:

        df[col] = pd.to_numeric(
            df[col],
            errors="coerce"
        )

    # --------------------------------------------------------
    # Check TOTAL
    # --------------------------------------------------------

    if "TOTAL" not in df.columns:

        raise ValueError(
            "TOTAL column not found in {}".format(
                filepath
            )
        )

    # Remove rows without TOTAL
    df = df.dropna(
        subset=["TOTAL"]
    ).copy()

    if len(df) == 0:

        raise ValueError(
            "No valid TOTAL values in {}".format(
                filepath
            )
        )

    return df


# ============================================================
# 3. ANALYZE EACH REPLICATE
# ============================================================

system_results = {}

replicate_results = {}


print("")
print("=" * 75)
print("MM/PBSA FRAME-BY-FRAME ANALYSIS")
print("=" * 75)
print("")


for system_name, file_list in systems.items():

    print("")
    print("-" * 75)
    print(system_name)
    print("-" * 75)

    replicate_values = []

    for rep_number, filepath in enumerate(
        file_list,
        start=1
    ):

        print(
            "Reading replicate {}: {}".format(
                rep_number,
                filepath
            )
        )

        df = load_mmpbsa_csv(filepath)

        # ----------------------------------------------------
        # Mean TOTAL across frames
        # ----------------------------------------------------

        mean_total = df["TOTAL"].mean()

        sd_total = df["TOTAL"].std(
            ddof=1
        )

        sem_total = stats.sem(
            df["TOTAL"]
        )

        n_frames = len(
            df["TOTAL"]
        )

        replicate_values.append(
            mean_total
        )

        print(
            "  Frames = {}".format(
                n_frames
            )
        )

        print(
            "  Mean TOTAL = {:.4f} kcal/mol".format(
                mean_total
            )
        )

        print(
            "  SD = {:.4f} kcal/mol".format(
                sd_total
            )
        )

        print(
            "  SEM = {:.4f} kcal/mol".format(
                sem_total
            )
        )

    # --------------------------------------------------------
    # Save replicate values
    # --------------------------------------------------------

    replicate_results[
        system_name
    ] = replicate_values

    # --------------------------------------------------------
    # Statistics across 3 independent replicas
    # --------------------------------------------------------

    mean_replicates = np.mean(
        replicate_values
    )

    sd_replicates = np.std(
        replicate_values,
        ddof=1
    )

    sem_replicates = stats.sem(
        replicate_values
    )

    system_results[
        system_name
    ] = {

        "Mean": mean_replicates,

        "SD": sd_replicates,

        "SEM": sem_replicates,

        "N_replicates": len(
            replicate_values
        )
    }

    print("")
    print(
        "Across replicas:"
    )

    print(
        "  Mean = {:.4f} kcal/mol".format(
            mean_replicates
        )
    )

    print(
        "  SD = {:.4f} kcal/mol".format(
            sd_replicates
        )
    )

    print(
        "  SEM = {:.4f} kcal/mol".format(
            sem_replicates
        )
    )


# ============================================================
# 4. GLOBAL SUMMARY
# ============================================================

summary_df = pd.DataFrame(
    system_results
).T


print("")
print("")
print("=" * 75)
print("GLOBAL MM/PBSA SUMMARY")
print("=" * 75)
print("")


print(
    summary_df[
        [
            "Mean",
            "SD",
            "SEM",
            "N_replicates"
        ]
    ].round(4)
)


summary_df.to_csv(
    "MM_PBSA_global_summary.csv"
)


# ============================================================
# 5. SAVE INDIVIDUAL REPLICATE VALUES
# ============================================================

replicate_rows = []


for system_name, values in replicate_results.items():

    for rep_number, value in enumerate(
        values,
        start=1
    ):

        replicate_rows.append({

            "System": system_name,

            "Replica": rep_number,

            "Mean_TOTAL_kcal_mol": value
        })


replicate_df = pd.DataFrame(
    replicate_rows
)


replicate_df.to_csv(
    "MM_PBSA_replicate_values.csv",
    index=False
)


# ============================================================
# 6. ETHANOL EFFECT
#
# DeltaDeltaG = +EtOH - -EtOH
# ============================================================

drugs = [
    "Mitragynine",
    "7OH",
    "TAN"
]


ddg_results = {}


for drug in drugs:

    no_etoh = np.array(
        replicate_results[
            "{}_-EtOH".format(drug)
        ],
        dtype=float
    )

    plus_etoh = np.array(
        replicate_results[
            "{}_+EtOH".format(drug)
        ],
        dtype=float
    )

    # --------------------------------------------------------
    # DeltaDeltaG
    # --------------------------------------------------------

    delta_delta_G = (
        np.mean(plus_etoh)
        - np.mean(no_etoh)
    )

    # --------------------------------------------------------
    # Paired t-test
    # --------------------------------------------------------

    t_stat, p_value = stats.ttest_rel(
        plus_etoh,
        no_etoh
    )

    ddg_results[
        drug
    ] = {

        "Minus_EtOH_Mean":
            np.mean(no_etoh),

        "Minus_EtOH_SD":
            np.std(
                no_etoh,
                ddof=1
            ),

        "Minus_EtOH_SEM":
            stats.sem(no_etoh),

        "Plus_EtOH_Mean":
            np.mean(plus_etoh),

        "Plus_EtOH_SD":
            np.std(
                plus_etoh,
                ddof=1
            ),

        "Plus_EtOH_SEM":
            stats.sem(plus_etoh),

        "DeltaDeltaG":
            delta_delta_G,

        "t_statistic":
            t_stat,

        "p_value":
            p_value
    }


# ============================================================
# 7. HOLM MULTIPLE-COMPARISON CORRECTION
# ============================================================

pvalues = [
    ddg_results[
        drug
    ]["p_value"]
    for drug in drugs
]


order = np.argsort(
    pvalues
)


adjusted = np.empty(
    len(pvalues)
)


for rank, index in enumerate(order):

    adjusted[index] = (
        pvalues[index]
        * (len(pvalues) - rank)
    )


# Force monotonicity
for i in range(1, len(order)):

    if adjusted[
        order[i]
    ] < adjusted[
        order[i - 1]
    ]:

        adjusted[
            order[i]
        ] = adjusted[
            order[i - 1]
        ]


adjusted = np.minimum(
    adjusted,
    1.0
)


for i, drug in enumerate(drugs):

    ddg_results[
        drug
    ]["p_adjusted_Holm"] = adjusted[i]

    ddg_results[
        drug
    ]["significant_Holm"] = (
        adjusted[i] < 0.05
    )


# ============================================================
# 8. ETHANOL SUMMARY TABLE
# ============================================================

ddg_df = pd.DataFrame(
    ddg_results
).T


print("")
print("")
print("=" * 75)
print("ETHANOL EFFECT")
print("=" * 75)
print("")


print(
    ddg_df.round(4)
)


ddg_df.to_csv(
    "MM_PBSA_ethanol_effect.csv"
)


# ============================================================
# 9. EASY-TO-READ ETHANOL RESULTS
# ============================================================

print("")
print("=" * 75)
print("ETHANOL EFFECT SUMMARY")
print("=" * 75)
print("")


for drug in drugs:

    result = ddg_results[
        drug
    ]

    print(
        "{}".format(drug)
    )

    print(
        "  -EtOH: {:.3f} +/- {:.3f} kcal/mol".format(
            result["Minus_EtOH_Mean"],
            result["Minus_EtOH_SD"]
        )
    )

    print(
        "  +EtOH: {:.3f} +/- {:.3f} kcal/mol".format(
            result["Plus_EtOH_Mean"],
            result["Plus_EtOH_SD"]
        )
    )

    print(
        "  DeltaDeltaG: {:.3f} kcal/mol".format(
            result["DeltaDeltaG"]
        )
    )

    print(
        "  Paired t-test p: {:.4f}".format(
            result["p_value"]
        )
    )

    print(
        "  Holm-adjusted p: {:.4f}".format(
            result["p_adjusted_Holm"]
        )
    )

    if result["significant_Holm"]:

        print(
            "  Statistical result: significant"
        )

    else:

        print(
            "  Statistical result: not significant"
        )

    print("")


# ============================================================
# 10. PLOT
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6),
    dpi=300
)


x = np.arange(
    len(drugs)
)

width = 0.35


# ------------------------------------------------------------
# Means
# ------------------------------------------------------------

means_no_etoh = [
    summary_df.loc[
        "{}_-EtOH".format(drug),
        "Mean"
    ]
    for drug in drugs
]


means_plus_etoh = [
    summary_df.loc[
        "{}_+EtOH".format(drug),
        "Mean"
    ]
    for drug in drugs
]


# ------------------------------------------------------------
# SEM
# ------------------------------------------------------------

sems_no_etoh = [
    summary_df.loc[
        "{}_-EtOH".format(drug),
        "SEM"
    ]
    for drug in drugs
]


sems_plus_etoh = [
    summary_df.loc[
        "{}_+EtOH".format(drug),
        "SEM"
    ]
    for drug in drugs
]


# ------------------------------------------------------------
# Bars
# ------------------------------------------------------------

ax.bar(
    x - width / 2,
    means_no_etoh,
    width,
    yerr=sems_no_etoh,
    capsize=5,
    label="- Ethanol",
    color="#2b5c8f",
    edgecolor="black",
    linewidth=0.8
)


ax.bar(
    x + width / 2,
    means_plus_etoh,
    width,
    yerr=sems_plus_etoh,
    capsize=5,
    label="+ Ethanol",
    color="#e08214",
    edgecolor="black",
    linewidth=0.8
)


# ============================================================
# 11. INDIVIDUAL REPLICATES
# ============================================================

rng = np.random.default_rng(
    42
)


for i, drug in enumerate(drugs):

    no_etoh = replicate_results[
        "{}_-EtOH".format(drug)
    ]

    plus_etoh = replicate_results[
        "{}_+EtOH".format(drug)
    ]

    jitter_no = rng.normal(
        0,
        0.025,
        len(no_etoh)
    )

    jitter_plus = rng.normal(
        0,
        0.025,
        len(plus_etoh)
    )

    ax.scatter(
        np.full(
            len(no_etoh),
            x[i] - width / 2
        ) + jitter_no,

        no_etoh,

        color="black",

        s=35,

        zorder=5
    )

    ax.scatter(
        np.full(
            len(plus_etoh),
            x[i] + width / 2
        ) + jitter_plus,

        plus_etoh,

        color="black",

        s=35,

        zorder=5
    )


# ============================================================
# 12. PLOT FORMATTING
# ============================================================

ax.set_ylabel(
    "MM/PBSA binding energy (kcal/mol)",
    fontsize=12,
    fontweight="bold"
)


ax.set_xticks(
    x
)


ax.set_xticklabels(
    [
        "Mitragynine",
        "7-OH-mitragynine",
        "TAN-67"
    ],
    fontsize=11,
    fontweight="bold"
)


ax.axhline(
    0,
    color="grey",
    linestyle="--",
    linewidth=0.8
)


ax.legend(
    frameon=True,
    facecolor="white",
    edgecolor="none"
)


ax.set_title(
    "Effect of Ethanol on Ligand-dOR MM/PBSA Energy",
    fontsize=13,
    fontweight="bold"
)


ax.spines[
    "top"
].set_visible(False)

ax.spines[
    "right"
].set_visible(False)


plt.tight_layout()


plt.savefig(
    "mmpbsa_ethanol_comparison.png",
    dpi=300,
    bbox_inches="tight"
)


plt.close()


# ============================================================
# 13. FINAL OUTPUT
# ============================================================

print("")
print("=" * 75)
print("ANALYSIS COMPLETE")
print("=" * 75)
print("")
print("Generated files:")
print("")
print("1. MM_PBSA_global_summary.csv")
print("2. MM_PBSA_replicate_values.csv")
print("3. MM_PBSA_ethanol_effect.csv")
print("4. mmpbsa_ethanol_comparison.png")
print("")