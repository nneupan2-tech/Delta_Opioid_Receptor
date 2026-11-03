Markdown
# Delta Opioid Receptor (DOR) Dynamics with Mitragynine, 7-Hydroxymitragynine and (−)-TAN-67 under Ethanol Co-solvent Exposure

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2025.2-blue.svg)](https://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
<!-- Uncomment post-publication / Zenodo upload:
[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Paper](https://img.shields.io/badge/JCIM-10.1021%2Facs.jcim.XXXXXXX-blue)](https://doi.org/10.1021/acs.jcim.XXXXXXX)
-->

This repository contains input structures, force-field parameters, execution protocols, analysis pipelines, and Markov State Modeling (MSM) workflows accompanying:

> **[TODO: Full Paper Title]**  
> *[TODO: Authors: e.g., Jane Doe, John Smith, ...]*  
> *Journal of Chemical Information and Modeling* (Submitted, 2026). DOI: [TODO: 10.1021/acs.jcim.XXXXXXX]  
> **Preprint:** [TODO: ChemRxiv DOI / link]

---

## Table of Contents

- [Overview](#overview)
- [System Architecture & Building](#system-architecture--building)
  - [Receptor Complex & Lipids](#receptor-complex--lipids)
  - [Ligand Preparation & Molecular Docking](#ligand-preparation--molecular-docking)
  - [Ethanol Co-solvent Insertion](#ethanol-co-solvent-insertion)
- [Simulation Matrix](#simulation-matrix)
- [Molecular Dynamics Protocol](#molecular-dynamics-protocol)
- [Post-Simulation Analysis Pipelines](#post-simulation-analysis-pipelines)
  - [Workflow Mapping](#workflow-mapping)
  - [Trajectory Preprocessing](#trajectory-preprocessing)
  - [Conformational & Activation Metrics](#conformational--activation-metrics)
  - [MM/PBSA Free Energy & Decomposition](#mmpbsa-free-energy--decomposition)
  - [TICA & Markov State Modeling (MSM)](#tica--markov-state-modeling-msm)
  - [Statistical Analysis & Visualization](#statistical-analysis--visualization)
- [Repository Structure](#repository-structure)
- [Data Availability & Reproduction](#data-availability--reproduction)
- [Citation & License](#citation--license)

---

## Overview

Ethanol exposure significantly alters G protein-coupled receptor (GPCR) dynamics and allosteric modulation. This project provides a comparative, all-atom molecular dynamics (MD) study of the active human **delta opioid receptor ($\delta$OR / DOR)** bound to three distinct agonists:

* **(−)-Mitragynine (MIT):** Primary alkaloid from *Mitragyna speciosa* (kratom).
* **7-Hydroxymitragynine (7OH):** Potent oxidized active metabolite of mitragynine.
* **(−)-TAN-67 (TAN):** Canonical, non-peptidic, highly selective DOR agonist.
* **Apo (Ligand-free):** Structural and baseline dynamic control.

Each condition was simulated in a physiological, heterogeneous lipid bilayer with and without **0.5% v/v ethanol co-solvent** across **3 independent replicates** for **1 µs each**, yielding **24 µs of cumulative all-atom sampling**.

---

## System Architecture & Building

### Receptor Complex & Lipids
* **Receptor Model:** Active-state cryo-EM structure of the human $\delta$-opioid receptor complexed with heterotrimeric $\text{G}_{i1}$ protein (PDB ID: [8Y45](https://www.rcsb.org/structure/8Y45); Cheng et al., 2024).
* **Included Subunits:** Chains A (DOR), B ($\text{G}\beta_1$), C ($\text{G}\gamma_2$), and E ($\text{G}\alpha_{i1}$).
* **Inactive Reference:** Inactive crystal structure of human DOR (PDB ID: [4N6H](https://www.rcsb.org/structure/4N6H); Fenalti et al., 2014).
* **Membrane Embedding:** Spatial orientation aligned using coordinates from the Orientations of Proteins in Membranes (OPM) database. Built via **CHARMM-GUI Bilayer Builder**.
* **Lipid Composition:** Asymmetric, physiologically relevant bilayer consisting of 450 total lipids (225 per leaflet):
  * **POPC** (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine): 102 molecules/leaflet
  * **POPE** (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine): 103 molecules/leaflet
  * **Cholesterol**: 20 molecules/leaflet
* **Solvation & Ionics:** Enclosed in a tetragonal periodic box ($\approx 115 \times 115 \times 165\text{ \AA}$) with TIP3P water molecules. Neutralized at an ionic strength of 0.15 M KCl using distance-based ion placement (153–156 $\text{K}^+$, 155–158 $\text{Cl}^-$).

### Ligand Preparation & Molecular Docking
* **2D Building & Geometry Optimization:** 2D topologies were constructed in ChemDraw Prime v22.2.0. Ligand energy minimization was carried out in MOE 2025 (Chemical Computing Group) to a gradient termination threshold of $0.1\text{ RMS kcal/mol/\AA}^2$.
* **Protonation:** The tertiary amine nitrogen was protonated at physiological pH 7.4; this charged state was conserved throughout parameterization and MD.
* **Docking Procedure:** MOE rigid-receptor docking targeted into the orthosteric biased-agonist (ADL5859) cavity of PDB 8Y45 using the triangle matcher placement method.
* **Scoring & Selection:** Initial screening with the London dG scoring function; refined via Generalized Born Volume Integral / Implicit Solvent (GBVI/WSA dG). Top poses exhibiting the canonical salt bridge between the protonated tertiary amine and **$\text{D}128^{3.32}$** were selected for MD system assembly.

| Ligand | Short ID | Residue Code | Protonation State | Parameterization |
| :--- | :--- | :--- | :--- | :--- |
| **(−)-Mitragynine** | MIT | `MIT` | $+1$ (tertiary amine protonated) | CHARMM General Force Field (CGenFF) |
| **7-Hydroxymitragynine** | 7OH | `7OH` | $+1$ (tertiary amine protonated) | CHARMM General Force Field (CGenFF) |
| **(−)-TAN-67** | TAN | `TAN` | $+1$ (tertiary amine protonated) | CHARMM General Force Field (CGenFF) |

### Ethanol Co-solvent Insertion
* **Concentration:** 0.5% v/v ($\approx 0.16\text{ mol\%}$), an experimentally and pharmacologically grounded concentration that avoids membrane denaturation while enabling microsecond kinetic observations.
* **Insertion:** Using `gmx insert-molecules`, **269 ethanol molecules** were randomly placed into the bulk aqueous phase of pre-equilibrated systems, replacing water molecules while strictly excluding the protein and bilayer core.

---

## Simulation Matrix

All simulations were run in three independent replicates (`r1`, `r2`, `r3`) using distinct pseudorandom initial velocity seeds (`gen-seed = -1`).

| System ID | Ligand | 0.5% v/v EtOH | Replicas $\times$ Length | Cumulative Time |
| :--- | :--- | :---: | :--- | :--- |
| `DOR_APO` | Apo (None) | No | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_MIT` | (−)-Mitragynine | No | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_7OH` | 7-Hydroxymitragynine | No | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_TAN` | (−)-TAN-67 | No | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_APO_EtOH` | Apo (None) | **Yes** | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_MIT_EtOH` | (−)-Mitragynine | **Yes** | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_7OH_EtOH` | 7-Hydroxymitragynine | **Yes** | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| `DOR_TAN_EtOH` | (−)-TAN-67 | **Yes** | $3 \times 1.0\text{ }\mu\text{s}$ | $3.0\text{ }\mu\text{s}$ |
| **Total** | | | **24 trajectories** | **$24.0\text{ }\mu\text{s}$** |

---

## Molecular Dynamics Protocol

All simulations were executed with **GROMACS 2025.2** compiled with CUDA 12.4 acceleration.

### Force Field & Integrator
* **Force Field:** CHARMM36m (proteins, lipids, ions), modified TIP3P (water), CGenFF v4.6 (ligands/ethanol).
* **Integrator:** Leap-frog algorithm with a 2 fs numerical integration time step.
* **Constraints:** All covalent bonds involving hydrogen atoms constrained via the LINCS algorithm.
* **Short-Range Interactions:** Non-bonded pair lists evaluated with the Verlet scheme. Van der Waals forces calculated using a smooth force-switch modifier from $1.0\text{ nm}$ to $1.2\text{ nm}$.
* **Long-Range Electrostatics:** Particle Mesh Ewald (PME) with a real-space grid cutoff of $1.2\text{ nm}$.

### Equilibration & Production Stages
1. **Minimization (`step6.0_minimization.mdp`):** 5,000 steps of steepest descent ($F_{\text{max}} < 1000\text{ kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$). Positional restraints ($1000\text{ kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$) applied to protein, ligand, and lipid heavy atoms.
2. **Stepwise Equilibration (`step6.1` – `step6.6_equilibration.mdp`):** Six-step progressive relaxation in NVT and NPT ensembles incrementally scaling down force constants on protein backbones, sidechains, lipid headgroups, and ligand atoms.
3. **Production (`step7_production.mdp`):**
   * **Ensemble:** Isothermal-isobaric ($NPT$), unrestrained for $1.0\text{ }\mu\text{s}$ per trajectory.
   * **Temperature Coupling:** 303.15 K using the Bussi-Donadio-Parrinello velocity-rescaling (`v-rescale`) thermostat ($\tau_t = 1.0\text{ ps}$). Independent coupling groups: `PROT`, `MEMB`, and `SOL_ION` (includes water, ions, ethanol, and ligand).
   * **Pressure Coupling:** 1.0 bar via semi-isotropic `c-rescale` barostat ($\tau_p = 5.0\text{ ps}$, compressibility $\beta = 4.5 \times 10^{-5}\text{ bar}^{-1}$).
   * **Output Frequency:** Trajectory frames recorded every $100\text{ ps}$ (10,000 frames/replica).

---

## Post-Simulation Analysis Pipelines

### Workflow Mapping

| Analysis | Script / Tool | Inputs & Window | Key Outputs |
| :--- | :--- | :--- | :--- |
| **PBC Correction & Fitting** | `analysis/preprocessing.sh` (`gmx trjconv`) | Raw `.xtc`, `topol.tpr` (full trajectory) | `traj_final_fixed.xtc` |
| **Global Stability & Dynamics** | `analysis/run_conformational.py` (`gmx rms`, `rmsfall`, `gyrate`) | Backbone / heavy atoms, 0–1000 ns | RMSD, RMSF, $R_g$ time series |
| **Intermolecular H-Bonds** | `analysis/hbond.sh` (`gmx hbond`), `analysis/hbond_occupancy.py` | Protein–ligand interface | `hbnum.xvg`, `hbmap_plot.png`, occupancy tables |
| **Receptor Activation Metrics** | `analysis/measure_d253_displacement.py`, `analysis/run_conformational.py` | TM3–TM6 $C_\alpha$, D253/W274, $\text{G}\alpha_{i1}$ $\alpha$5-helix | Inter-helix distance, microswitch shift, G-protein RMSD |
| **Binding Free Energy (MM/PBSA)** | `analysis/mmpbsa_run.sh` (`gmx_MMPBSA` v1.6+ / AmberTools) | Final 250 ns (frames 750–1000), 2-frame stride | `FINAL_RESULTS_GIBBS_MMPBSA.dat`, per-residue decomposition |
| **Kinetics & Landscapes (MSM)** | `analysis/run_msm_deeptime.py` (`deeptime`) | 3D feature vectors ($\tau_{\text{TICA}}=1.0\text{ ns}$, $\tau_{\text{MSM}}=1.5\text{ ns}$) | 2D TICA projections, PCCA+ macrostates, CK-test plots |

### Trajectory Preprocessing
Prior to metric evaluation, all raw trajectories were centered, corrected for periodic boundary conditions (PBC), and aligned against the starting structure:

```bash
# 1. Cluster whole molecules and remove periodic jumps
gmx trjconv -s topol.tpr -f traj_comp.xtc -o traj_noPBC.xtc -pbc nojump

# 2. Center receptor complex in the box
gmx trjconv -s topol.tpr -f traj_noPBC.xtc -o traj_centered.xtc -pbc mol -center

# 3. Fit rotational/translational motions to the C-alpha atoms
gmx trjconv -s topol.tpr -f traj_centered.xtc -o traj_final_fixed.xtc -fit rot+trans
Conformational & Activation Metrics
Calculated via analysis/run_conformational.py and standard GROMACS analysis suites:

Global Stability: RMSD (backbone and ligand), RMSF, and radius of gyration (
R 
g
​
 ).

Hydrogen Bond Networks: Intermolecular ligand-protein hydrogen bonding frequencies and occupancies (gmx hbond).

Receptor Activation Metrics:

TM6 Outward Tilting: Displacement of TM6 relative to the inactive crystal structure (PDB ID: 4N6H).

TM3–TM6 Inter-Helix Distance: Distance between 
C 
α
​
  atoms of TM3 and TM6 to quantify cytoplasmic intracellular opening.

Micro-switch Motion: Positional displacement of residue D253 / W274 across the trajectories relative to inactive 4N6H.

G-Protein Coupling Stability: Root-mean-square deviation of the C-terminal 
α5-helix of the 
Gα 
i1
​
  subunit.

MM/PBSA Free Energy & Decomposition
Evaluated using the gmx_MMPBSA v1.6+ suite interfaced with AmberTools sander:

Sampled Window: The final 
250 ns of production MD (frames 750–1000) sampled at a 2-frame interval.

Implicit Solvent / Membrane PB Model:

ΔG 
bind
​
 =ΔE 
MM
​
 +ΔG 
solv
​
 −TΔS≈ΔE 
vdW
​
 +ΔE 
elec
​
 +ΔG 
PB
​
 +ΔG 
NP
​
 
Dielectric Constants: Solute 
ε 
in
​
 =2.0, Aqueous Phase 
ε 
out
​
 =80.0, Membrane Bilayer 
ε 
mem
​
 =2.0.

Membrane Parameters: Implicit membrane thickness 
=35.0  
A
˚
 , ionic strength 
=0.15 M.

Entropy: Configurational entropy terms (
TΔS) were omitted due to the high structural homology among compared ligands, approximating relative affinities via binding enthalpy.

Per-Residue Decomposition: Backbone and sidechain per-residue decomposition within a 
5.0  
A
˚
  sphere of the binding pocket (energetic threshold 
∣ΔG∣≥0.5 kcal/mol).

TICA & Markov State Modeling (MSM)
Kinetic transitions and free-energy landscapes were constructed using the deeptime library:

Feature Construction: Concatenated structural order parameters into a 3D feature vector 
x(t):

TM3–TM6 
C 
α
​
  inter-helix distance.

Microswitch residue displacement relative to 4N6H.

G-protein C-terminal 
α5-helix backbone RMSD.

Time-Lagged Independent Component Analysis (TICA):

Applied globally across concatenated replicates with a lag time 
τ 
TICA
​
 =1.0 ns (10 frames) to resolve slow relaxation modes:

C(0)v 
i
​
 =λ 
i
​
 C(τ)v 
i
​
 
Clustering & Discretization: The 2D TICA projection was partitioned into 
k=100 discrete microstates using regular 
k-means clustering.

Transition Matrix & Macrostates:

Maximum Likelihood MSMs were estimated at an MSM lag time of 
τ 
MSM
​
 =1.5 ns (15 frames) enforcing detailed balance.

Perron Cluster Cluster Analysis (PCCA+): Grouped microstates into 3 primary macrostates:

M0: Uncoupled / Inactive-like

M1: Pre-coupled intermediate

M2: Fully coupled active state

Validation:

Implied Timescales (
t 
i
​
 ): Computed across varying lag times 
τ to confirm kinetic convergence and memoryless Markovianity.

Chapman–Kolmogorov Test: Evaluated to ensure trajectory dynamics match model propagation:

T(kτ)≈T(τ) 
k
 
Statistical Analysis & Visualization
Visualization Tools: Three-dimensional binding modes, membrane orientations, and electrostatic surface potentials were visualized using PyMOL v3.1.8 and Discovery Studio Visualizer 4.5.

Statistical Rigor: Statistical comparisons were conducted across independent simulation replicates (
N=3 replicates per experimental condition). Pairwise differences between conditions were evaluated using Welch's independent samples unequal-variance 
t-test.

Repository Structure
Plaintext
.
├── LICENSE                          # MIT open-source code license
├── CITATION.cff                     # Citation metadata
├── environment.yml                  # Conda environment configuration
├── data/                            # System setups and input definitions
│   ├── DOR_APO/
│   │   ├── shared/                  # Common topology, parameters, and inputs
│   │   │   ├── topol.top            # GROMACS topology
│   │   │   ├── toppar/              # CHARMM36m & CGenFF parameters
│   │   │   ├── inputs/              # step6.0_minimization to step7_production .mdp
│   │   │   └── run_gromacs.sh       # Execution wrapper
│   │   └── replicas/
│   │       ├── r1/                  # Production logs, seeds, and run outputs
│   │       ├── r2/
│   │       └── r3/
│   ├── DOR_MIT/
│   ├── DOR_7OH/
│   ├── DOR_TAN/
│   ├── DOR_APO_EtOH/
│   ├── DOR_MIT_EtOH/
│   ├── DOR_7OH_EtOH/
│   └── DOR_TAN_EtOH/
├── analysis/                        # Production analysis workflows
│   ├── preprocessing.sh             # PBC removal, centering, and alignment
│   ├── run_conformational.py        # RMSD, RMSF, Rg, and G-protein displacement
│   ├── measure_d253_displacement.py # Microswitch D253 displacement
│   ├── mmpbsa_run.sh                # gmx_MMPBSA execution script
│   ├── mmpbsa.in                    # MMPBSA input parameters
│   └── run_msm_deeptime.py          # TICA, MSM, PCCA+, and CK-test script
├── results/                         # Processed time-series and output plots
│   ├── mmpbsa/                      # Free energy tables and decomposition CSVs
│   ├── activation_metrics/          # Inter-helix distance and RMSD plots
│   └── msm/                         # Implied timescales, free energy contours
└── trajectories/
    └── MANIFEST.md                  # Zenodo DOIs, SHA-256 checksums, and file sizes
Data Availability & Reproduction
1. Environment Setup
Install dependencies via Conda:

Bash
conda env create -f environment.yml
conda activate dor-md
2. Large File Retrieval
Because raw trajectories exceed repository storage limits, .xtc, .tpr, and coordinate files are archived on Zenodo:

Bash
# Retrieve production trajectories and binary files
bash scripts/download_zenodo_data.sh

# Verify download integrity
sha256sum -c trajectories/SHA256SUMS.txt
3. Running Analysis Workflows
To reproduce conformational and activation metrics across replicas:

Bash
# Execute conformational and activation metrics for DOR_7OH
python analysis/run_conformational.py --system DOR_7OH --replicas r1 r2 r3

# Execute MM/PBSA binding free energy calculation
bash analysis/mmpbsa_run.sh DOR_7OH

# Run TICA and Markov State Modeling
python analysis/run_msm_deeptime.py --lag_tica 10 --lag_msm 15
Citation & License
Citation
If you utilize this workflow, topologies, or simulation trajectories in your research, please cite:

Code snippet
@article{dor_ethanol_dynamics_2026,
  author    = {[TODO: Authors]},
  title     = {[TODO: Full paper title]},
  journal   = {J. Chem. Inf. Model.},
  year      = {2026},
  doi       = {[TODO: 10.1021/acs.jcim.XXXXXXX]}
}
License
Code and Workflows (analysis/, data/**/*.sh, scripts/): Released under the MIT License.

Data, Topologies, and Processed Trajectories (data/, results/, and Zenodo records): Distributed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

Contact
Primary Investigator: [TODO: Name] ([TODO: Email])

Laboratory / Institution: [TODO: Research Group, University / Institute]

ORCID: [TODO: XXXX-XXXX-XXXX-XXXX]
