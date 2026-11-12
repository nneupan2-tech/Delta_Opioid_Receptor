# Delta Opioid Receptor (DOR) Dynamics with Mitragynine, 7-Hydroxymitragynine and (−)-TAN-67 under Ethanol Co-solvent Exposure

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![GROMACS](https://img.shields.io/badge/GROMACS-2025.2-blue.svg)](https://www.gromacs.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)

This repository contains input structures, force-field parameters, execution protocols, analysis pipelines, and Markov State Modeling (MSM) workflows accompanying:

> **Delta Opioid Receptor (DOR) Dynamics with Mitragynine, 7-Hydroxymitragynine and (−)-TAN-67 under Ethanol Co-solvent Exposure**  
> *Manuscript in preparation (2026).*  
> *(Preprint and data repository links will be updated upon deposit).*

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
  - [MM/PBSA Free Energy & Decomposition Analysis](#mmpbsa-free-energy--decomposition-analysis)
  - [TICA & Markov State Modeling (MSM)](#tica--markov-state-modeling-msm)
  - [Statistical Analysis & Visualization](#statistical-analysis--visualization)
- [Repository Structure](#repository-structure)
- [Data Availability & Reproduction](#data-availability--reproduction)
  - [1. Environment Setup](#1-environment-setup)
  - [2. Large File Retrieval](#2-large-file-retrieval)
  - [3. Running Analysis Workflows](#3-running-analysis-workflows)
- [Citation & License](#citation--license)

---

## Overview

Ethanol exposure significantly modulates G protein-coupled receptor (GPCR) conformations, altering basal signaling and allosteric communication. This repository contains the complete computational workflow for an all-atom molecular dynamics (MD) study probing the human **delta opioid receptor ($\delta$OR / DOR)** bound to three structurally distinct agonists alongside an unliganded control:

* **(−)-Mitragynine (MIT):** Principal psychoactive alkaloid derived from *Mitragyna speciosa* (kratom).
* **7-Hydroxymitragynine (7OH):** Highly potent, oxidized active metabolite of mitragynine.
* **(−)-TAN-67 (TAN):** Non-peptidic, highly selective standard $\delta$OR agonist.
* **Apo (Ligand-free):** Unliganded control used to determine baseline receptor plasticity.

Each condition was simulated in a complex, heterogeneous lipid bilayer in both the absence and presence of **$0.5\%\text{ v/v}$ ethanol co-solvent** across **3 independent replicates** for **$1.0\text{ }\mu\text{s}$ each**, totaling **$24.0\text{ }\mu\text{s}$ of cumulative sampling**.

### Where to Find What

| Content | Location | Description |
| :--- | :--- | :--- |
| **Code & MDPs** | This repository | Simulation parameter files (`.mdp`), analysis wrappers, Python scripts |
| **Topologies** | `data/<system>/shared/` | System coordinates, `.top`, `.itp`, and parameter definitions |
| **Raw Trajectories** | Zenodo (10.5281/zenodo.23955575.) | Production `.xtc`, `.tpr`, `.gro`, and checkpoint files |
| **Manuscript** | TODO | Peer-reviewed publication and supplementary dataset |

---

## System Architecture & Building

### Receptor Complex & Lipids
* **Receptor Starting Model:** Active-state cryo-EM structure of human $\delta$OR complexed with the heterotrimeric $\text{G}_{\text{i}1}$ protein (PDB ID: [8Y45](https://www.rcsb.org/structure/8Y45); Cheng et al., 2024).
* **Included Subunits:** Chains A (DOR), B ($\text{G}\beta_1$), C ($\text{G}\gamma_2$), and E ($\text{G}\alpha_{\text{i}1}$).
* **Inactive Reference:** High-resolution crystal structure of inactive human $\delta$OR (PDB ID: [4N6H](https://www.rcsb.org/structure/4N6H); Fenalti et al., 2014).
* **Membrane Embedding:** Spatially positioned within the bilayer using coordinates provided by the Orientations of Proteins in Membranes (OPM) database. Systems were assembled using the **CHARMM-GUI Bilayer Builder**.
* **Lipid Bilayer Composition:** Asymmetric bilayer consisting of 450 total lipids (225 molecules per leaflet):
  * **POPC** (1-palmitoyl-2-oleoyl-$sn$-glycero-3-phosphocholine): 102 molecules/leaflet
  * **POPE** (1-palmitoyl-2-oleoyl-$sn$-glycero-3-phosphoethanolamine): 103 molecules/leaflet
  * **Cholesterol**: 20 molecules/leaflet
* **Solvation & Ionics:** Enclosed within a tetragonal periodic boundary box ($\approx 115 \times 115 \times 165\text{ \AA}$) with TIP3P water molecules. Electrostatic neutralization and physiological ionic strength ($0.15\text{ M KCl}$) were achieved using distance-based ion placement (153–156 $\text{K}^+$, 155–158 $\text{Cl}^-$).

### Ligand Preparation & Molecular Docking
* **Ligand Topology & Optimization:** 2D structures of MIT, 7OH, and TAN were built in ChemDraw Prime v22.2.0. Energy minimization was performed in MOE 2025 (Chemical Computing Group) to an RMS gradient convergence threshold of $0.1\text{ kcal}\cdot\text{mol}^{-1}\cdot\text{\AA}^{-2}$.
* **Protonation State:** The aliphatic tertiary amine was protonated at physiological $\text{pH } 7.4$ ($+1$ net charge state), conserved across force-field parameter generation and subsequent dynamics.
* **Docking Protocol:** Rigid-receptor docking in MOE targeting the orthosteric biased-agonist (ADL5859) cavity of PDB 8Y45 via triangle matcher placement. Poses were scored using the London $\text{d}G$ scoring function and refined with the Generalized Born Volume Integral / Implicit Solvent (GBVI/WSA $\text{d}G$) free energy estimator.
* **Pose Filtering:** Docked poses displaying the canonical salt bridge interaction between the protonated ligand amine and the conserved aspartate **$\text{D}128^{3.32}$** were selected for molecular dynamics system assembly.

| Ligand | Short ID | Residue Code | Charge State | Force Field Parameterization |
| :--- | :--- | :--- | :---: | :--- |
| **(−)-Mitragynine** | MIT | `MIT` | $+1$ | CHARMM General Force Field (CGenFF v4.6) |
| **7-Hydroxymitragynine** | 7OH | `7OH` | $+1$ | CHARMM General Force Field (CGenFF v4.6) |
| **(−)-TAN-67** | TAN | `TAN` | $+1$ | CHARMM General Force Field (CGenFF v4.6) |

### Ethanol Co-solvent Insertion
* **Target Concentration:** $0.5\%\text{ v/v}$ ($\approx 0.16\text{ mol\%}$), approximating pharmacologically achievable concentrations while retaining sufficient cosolvent sampling within the microsecond regime.
* **Solvent Replacement:** Using `gmx insert-molecules`, 269 ethanol molecules were randomly inserted into the bulk aqueous phase of the equilibrated receptor-bilayer system, substituting water molecules while strictly excluding the protein interior and hydrophobic core of the bilayer.

---

| System ID | Ligand | 0.5% (v/v) EtOH | Replicas × Length | Cumulative Sampling |
| :--- | :--- | :---: | :---: | :---: |
| **DOR_APO** | None (Apo) | No | 3 × 1.0 μs | 3.0 μs |
| **DOR_MIT** | (−)-Mitragynine | No | 3 × 1.0 μs | 3.0 μs |
| **DOR_7OH** | 7-Hydroxymitragynine | No | 3 × 1.0 μs | 3.0 μs |
| **DOR_TAN** | (−)-TAN-67 | No | 3 × 1.0 μs | 3.0 μs |
| **DOR_APO_EtOH** | None (Apo) | Yes | 3 × 1.0 μs | 3.0 μs |
| **DOR_MIT_EtOH** | (−)-Mitragynine | Yes | 3 × 1.0 μs | 3.0 μs |
| **DOR_7OH_EtOH** | 7-Hydroxymitragynine | Yes | 3 × 1.0 μs | 3.0 μs |
| **DOR_TAN_EtOH** | (−)-TAN-67 | Yes | 3 × 1.0 μs | 3.0 μs |
| **Total** | | | **24 trajectories** | **24.0 μs** |

### Velocity Random Seeds

Initial atomic velocities for replicas r1, r2, and r3 were independently assigned from a Maxwell–Boltzmann distribution at 303.15 K using pseudorandom generation (gen-vel = yes, gen-seed = -1) in GROMACS to ensure statistically independent trajectories.

## Molecular Dynamics Protocol

All simulations were run on **GROMACS 2025.2** accelerated via CUDA 12.4.

### Force Field & Integrator Settings
* **Force Field:** CHARMM36m for protein, lipid bilayer, and aqueous ions; modified TIP3P for water; CGenFF v4.6 for ligands and ethanol.
* **Numerical Integrator:** Leap-frog algorithm with an integration time step of $\Delta t = 2\text{ fs}$.
* **Constraint Algorithm:** All covalent bonds containing hydrogen atoms were constrained via LINCS.
* **Short-Range Cutoffs:** Non-bonded pair lists generated using the Verlet scheme. Van der Waals forces smoothed via a force-switch function from $1.0\text{ nm}$ to $1.2\text{ nm}$.
* **Electrostatics:** Particle Mesh Ewald (PME) with a real-space cutoff of $1.2\text{ nm}$ and a Fourier grid spacing of $0.12\text{ nm}$.

### Staged Equilibration & Production
1. **Energy Minimization (`step6.0_minimization.mdp`):** 5,000 steps of steepest descent until $F_{\text{max}} < 1000\text{ kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$. Harmonic restraints ($1000\text{ kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$) applied to all non-hydrogen atoms of protein, ligand, and lipids.
2. **Equilibration Cascade (`step6.1` – `step6.6_equilibration.mdp`):** Six-step progressive relaxation in $NVT$ and $NPT$ ensembles systematically reducing positional restraints on protein backbone, sidechains, lipid headgroups, and ligand atoms.
3. **Production MD (`step7_production.mdp`):**
   * **Ensemble:** Isothermal-isobaric ($NPT$) without positional restraints for $1.0\text{ }\mu\text{s}$ per replica.
   * **Temperature Coupling:** Maintained at $303.15\text{ K}$ via the Bussi-Donadio-Parrinello velocity-rescaling (`v-rescale`) thermostat ($\tau_t = 1.0\text{ ps}$), coupled independently to `PROT`, `MEMB`, and `SOL_ION` (water, ions, ligand, and ethanol).
   * **Pressure Coupling:** Maintained at $1.0\text{ bar}$ via semi-isotropic `c-rescale` barostat ($\tau_p = 5.0\text{ ps}$, compressibility $\beta = 4.5 \times 10^{-5}\text{ bar}^{-1}$).
   * **Output Frequency:** Trajectory frames recorded every $100\text{ ps}$ (10,000 frames recorded per replicate).

---

## Post-Simulation Analysis Pipelines

### Workflow Mapping

| Analysis Module | Executable / Tool | Target Selection & Interval | Output Artifacts |
| :--- | :--- | :--- | :--- |
| **PBC Correction & Alignment** | `analysis/preprocessing.sh` (`gmx trjconv`) | Full system, whole trajectory | `traj_final_fixed.xtc` |
| **Global Conformational Dynamics** | `analysis/run_conformational.py` (`gmx rms`, `rmsfall`, `gyrate`) | Receptor $C_\alpha$, ligand heavy atoms ($0\text{--}1000\text{ ns}$) | RMSD, RMSF, $R_g$ time-series (`results/activation_metrics/`) |
| **H-Bond Analysis & Occupancy** | `analysis/hbond.sh` (`gmx hbond`), `analysis/hbond_occupancy.py` | Protein–ligand interface donors/acceptors | `hbnum.xvg`, `hbmap_plot.png`, occupancy tables |
| **Receptor Activation Metrics** | `analysis/measure_d253_displacement.py`, `analysis/run_conformational.py` | TM3–TM6 distance, D253/W274 shifts, $\text{G}\alpha_{\text{i}1}$ $\alpha$5-helix | Intracellular distance profiles, microswitch vectors, G-protein RMSD |
| **Binding Free Energy (MM/PBSA)** | `analysis/mmpbsa_run.sh` (`gmx_MMPBSA` v1.6+ / AmberTools) | Production frames 750–1000 (stride = 2 frames) | `FINAL_RESULTS_GIBBS_MMPBSA.dat`, per-residue decomposition tables |
| **Kinetic Landscapes (MSM)** | `analysis/run_msm_deeptime.py` (`deeptime`) | 3D feature vectors ($\tau_{\text{TICA}} = 1.0\text{ ns}$, $\tau_{\text{MSM}} = 1.5\text{ ns}$) | 2D TICA projections, PCCA+ macrostate kinetics, CK-test plots |

---

### Trajectory Preprocessing
Raw trajectory files were processed to remove periodic boundary jumps, center the macromolecular complex, and fit rotational/translational frames to the starting structure:

# 1. Cluster intact molecules across periodic boundaries
gmx trjconv -s topol.tpr -f traj_comp.xtc -o traj_noPBC.xtc -pbc nojump

# 2. Center receptor complex within the simulation box
gmx trjconv -s topol.tpr -f traj_noPBC.xtc -o traj_centered.xtc -pbc mol -center

# 3. Fit rotational and translational motions to receptor C-alpha atoms
gmx trjconv -s topol.tpr -f traj_centered.xtc -o traj_final_fixed.xtc -fit rot+trans

---

### Conformational & Activation Metrics
Calculated via `analysis/run_conformational.py` and standard GROMACS tools:
* **Global Structural Stability:** Root-mean-square deviation (RMSD) of receptor backbone and ligand heavy atoms, per-residue root-mean-square fluctuation (RMSF), and radius of gyration ($R_g$).
* **Hydrogen Bond Persistence:** Intermolecular ligand–receptor hydrogen bond counts and fractional occupancies computed using `gmx hbond`.
* **Receptor Activation & Coupling Metrics:**
  * **TM6 Outward Tilting:** Outward displacement of Transmembrane Helix 6 measured relative to inactive DOR (PDB ID: [4N6H](https://www.rcsb.org/structure/4N6H)).
  * **TM3–TM6 Inter-Helix Distance:** Distance between the $C_\alpha$ atoms of TM3 and TM6 to quantify cytoplasmic intracellular opening.
  * **Micro-switch Motion:** Positional displacement of residue **D253** (and adjacent switch residue **W274**) tracked against inactive 4N6H.
  * **G-Protein Coupling Stability:** Relative RMSD of the C-terminal $\alpha$5-helix of the $\text{G}\alpha_{\text{i}1}$ subunit.

---

### MM/PBSA Free Energy & Decomposition Analysis
Evaluated using `gmx_MMPBSA` (v1.6+) integrated with AmberTools `sander`:
* **Evaluation Window:** Sampled across the terminal $250\text{ ns}$ of production MD (frames 750 to 1000) at an interval of every two frames (125 frames/replica; 375 frames/system).
* **Implicit Solvent / Membrane PB Formulation:**
  The net binding free energy was partitioned into molecular mechanical and continuum solvation components:
  $$\Delta G_{\text{bind}} = \Delta E_{\text{MM}} + \Delta G_{\text{solv}} - T\Delta S \approx \Delta E_{\text{vdW}} + \Delta E_{\text{elec}} + \Delta G_{\text{PB}} + \Delta G_{\text{NP}}$$
  where:
  $$\Delta E_{\text{MM}} = \Delta E_{\text{internal}} + \Delta E_{\text{electrostatic}} + \Delta E_{\text{vdW}}$$
  $$\Delta G_{\text{solv}} = \Delta G_{\text{PB}} + \Delta G_{\text{NP}}$$
* **Dielectric & Membrane Parameters:**
  * Solute internal dielectric constant: $\varepsilon_{\text{in}} = 2.0$
  * External aqueous solvent dielectric: $\varepsilon_{\text{out}} = 80.0$
  * Hydrophobic membrane slab dielectric: $\varepsilon_{\text{mem}} = 2.0$
  * Implicit membrane thickness: $35.0\text{ \AA}$ (centered at the bilayer midplane)
  * Continuum ionic strength: $0.15\text{ M}$
* **Entropy Term:** Conformational entropy ($-T\Delta S$) was omitted due to the structural homology among the ligands, allowing relative affinities to be approximated by binding enthalpy.
* **Per-Residue Decomposition:** Backbone and sidechain electrostatic, van der Waals, and polar solvation contributions were computed for all protein residues within a $5.0\text{ \AA}$ sphere of the ligand binding pocket (significance threshold: $|\Delta G| \ge 0.5\text{ kcal}\cdot\text{mol}^{-1}$).

---

### TICA & Markov State Modeling (MSM)
Kinetic transitions and free-energy landscapes were constructed using the `deeptime` library:

1. **Feature Construction:** Concatenated structural order parameters into a 3D feature vector $\mathbf{x}(t)$:
   * TM3–TM6 $C_\alpha$ inter-helix distance.
   * Positional displacement of microswitch residue D253 relative to 4N6H.
   * Relative RMSD of the $\text{G}\alpha_{\text{i}1}$ C-terminal $\alpha$5-helix.
2. **Time-Lagged Independent Component Analysis (TICA):**
   * Applied globally across concatenated replicates with a lag time $\tau_{\text{TICA}} = 1.0\text{ ns}$ (10 frames) to resolve slow relaxation modes:
     $$C(0)\mathbf{v}_i = \lambda_i C(\tau)\mathbf{v}_i$$
     where $C(0)$ and $C(\tau)$ are the instantaneous and time-lagged covariance matrices.
3. **Clustering & Microstate Discretization:** The continuous 2D TICA space was discretized into $k = 100$ microstate clusters using $k$-means clustering.
4. **Transition Matrix Estimation:** Reversible maximum-likelihood transition matrices were estimated at an MSM lag time of $\tau_{\text{MSM}} = 1.5\text{ ns}$ (15 frames) enforcing detailed balance:
   $$\pi_i P_{ij}(\tau) = \pi_j P_{ji}(\tau)$$
5. **Macrostates (PCCA+):** Microstates were clustered into three functional macrostates using Perron Cluster Cluster Analysis:
   * **M0:** Uncoupled / Inactive-like state
   * **M1:** Pre-coupled intermediate state
   * **M2:** Fully coupled, active agonist-bound state
6. **Free Energy Landscapes & State Properties:**
   * **2D Potential of Mean Force:**
     $$F_i = -k_{\text{B}} T \ln \pi_i, \quad \text{normalized to } \min(F_i) = 0$$
   * **Thermodynamically Weighted Feature Vectors:**
     $$\bar{\mathbf{x}}_\alpha = \sum_{i} \chi_{\alpha i} \mathbf{c}_i$$
     where $\chi_{\alpha i}$ is the PCCA+ membership degree of microstate $i$ in macrostate $\alpha$, and $\mathbf{c}_i$ is the cluster center in feature space.
7. **Model Validation:**
   * **Implied Timescales ($t_i$):** Computed across lag times up to $\tau = 5.0\text{ ns}$ ($t_i = -\tau / \ln |\lambda_i|$) to verify timescale plateauing and memoryless Markovianity.
   * **Chapman–Kolmogorov Test:** Evaluated by propagating the estimated transition matrix against trajectory observations:
     $$P(k\tau) \approx P(\tau)^k$$

---

### Statistical Analysis & Visualization
* **3D Molecular Graphics:** Receptor-ligand interactions, membrane positioning, and electrostatic potential surfaces were generated using PyMOL v3.1.8 and Discovery Studio Visualizer 4.5.
* **Statistical Rigor:** Statistical comparisons were conducted across independent simulation replicates ($N = 3$ per condition). Pairwise differences between conditions (e.g., control vs. $+0.5\%$ ethanol) were evaluated via Welch's two-sample unequal-variance $t$-test ($p < 0.05$ considered statistically significant).

---

## Repository Structure

```text
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
```

---

## Data Availability & Reproduction

# Retrieve production trajectories and binary files
bash scripts/download_zenodo_data.sh

# Verify download integrity
sha256sum -c trajectories/SHA256SUMS.txt

# Execute conformational and activation metrics for DOR_7OH
python analysis/run_conformational.py --system DOR_7OH --replicas r1 r2 r3

# Execute MM/PBSA binding free energy calculation
bash analysis/mmpbsa_run.sh DOR_7OH

# Run TICA and Markov State Modeling
python analysis/msm.py --lag_tica 10 --lag_msm 15

---

## Citation & License

### Citation
If you utilize this workflow, topologies, or simulation trajectories in your research, please cite:
```bibtex
@article{dor_ethanol_dynamics_2026,
  author    = {[TODO: Authors]},
  title     = {[TODO: Full paper title]},
  journal   = {J. Chem. Inf. Model.},
  year      = {2026},
  doi       = {[TODO: 10.1021/acs.jcim.XXXXXXX]}
}
```

### License
* **Code and Workflows** (`analysis/`, `data/**/*.sh`, `scripts/`): Released under the [MIT License](LICENSE).
* **Data, Topologies, and Processed Trajectories** (`data/`, `results/`, and Zenodo records): Distributed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

### Contact
* **Primary Investigator:** Joydip Das 
* **Laboratory / Institution:** Dr Das Lab, University of Houston
* **ORCID:**  https://orcid.org/0000-0003-1540-5804
