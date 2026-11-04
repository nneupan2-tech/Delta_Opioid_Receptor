#!/usr/bin/env python3
"""
H-bond occupancy per residue from GROMACS gmx hbond outputs (-hbn .ndx, -hbm .xpm).

Expected files in DATA_DIR:
    hbond_<system>_r<N>.ndx
    hbmap_<system>_r<N>.xpm
    1000ns_<system>_r1.pdb     (atom order = atom numbering used in the ndx)

Outputs (in OUT_DIR):
    hbond_atom_level.csv        every H-bond, every replicate (atom names + occupancy)
    hbond_residue_level.csv     residue-pair occupancy per replicate
    hbond_summary.csv           mean +/- SD over replicates, all residues
    table_main_>=20pct.csv      residues with mean occupancy >= MAIN_CUTOFF in any system
    fig_main_heatmap.png        main-text figure (all systems side by side)
    supp/bar_<system>.png       one supplementary figure per system
"""
import os, re, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------- settings -----------------------------
DATA_DIR = "."
OUT_DIR = "hbond_results"
SYSTEMS = ["7oh", "mit", "tan", "7ohetoh", "mitetoh", "tanetoh"]
MAIN_CUTOFF = 20.0     # % : residues shown in main table
HEATMAP_CUTOFF = 10.0  # % : residues shown in main heatmap
SUPP_CUTOFF = 5.0      # % : residues shown in supplementary bar plots
REVERSE_XPM_ROWS = True  # GROMACS writes xpm rows top = last H-bond

PROTEIN = set("""ALA ARG ASN ASP CYS CYX CYM GLN GLU GLH GLY HIS HID HIE HIP HSD HSE HSP
ILE LEU LYS LYN MET PHE PRO SER THR TRP TYR VAL ASH NALA NGLY ACE NME""".split())

# ----------------------------- parsers -----------------------------
def read_pdb(path):
    """Return DataFrame indexed 1..N (by order of atoms, first model only)."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                rows.append((line[12:16].strip(), line[17:21].strip(),
                             line[21].strip(), int(line[22:26])))
            elif line.startswith("ENDMDL"):
                break
    df = pd.DataFrame(rows, columns=["atom", "resname", "chain", "resnum"])
    df.index = np.arange(1, len(df) + 1)
    return df

def read_ndx_triplets(path, n_rows):
    """Return array (n_rows,3) of 1-based atom numbers: donor, hydrogen, acceptor.
    The ndx also contains plain atom groups (e.g. [ Protein ]); we pick the group
    that has exactly 3*n_rows entries (the number of H-bonds in the xpm)."""
    groups, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):
                name = line.strip("[] ").strip()
                groups[name] = []
            elif line and name is not None:
                groups[name] += [int(x) for x in line.split()]
    cands = [k for k, v in groups.items() if len(v) == 3 * n_rows]
    if not cands:
        raise ValueError(f"{path}: no group with {3*n_rows} entries. Groups: "
                         f"{ {k: len(v) for k, v in groups.items()} }")
    key = next((k for k in cands if "hbond" in k.lower()), cands[0])
    if len(cands) > 1:
        print(f"  note: several candidate groups {cands}, using [{key}]")
    return np.array(groups[key]).reshape(-1, 3)

def read_xpm(path):
    """Return boolean array (n_hbonds, n_frames); True = H-bond present."""
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    q = [l for l in lines if l.lstrip().startswith('"')]
    w, h, nc, cpp = [int(x) for x in q[0].split('"')[1].split()]
    present_char = None
    for l in q[1:1 + nc]:
        s = l.split('"')[1]
        ch, rest = s[:cpp], s[cpp:]
        comment = l.split("/*")[-1].lower() if "/*" in l else ""
        if "present" in comment or ("yes" in comment) or ("hbond" in comment and "none" not in comment):
            present_char = ch
        elif present_char is None and "#ffffff" not in rest.lower() and "none" not in comment:
            present_char = ch
    data = [l.split('"')[1] for l in q[1 + nc:1 + nc + h]]
    mat = np.array([[c == present_char for c in row] for row in data])
    return mat[::-1] if REVERSE_XPM_ROWS else mat

# ----------------------------- analysis -----------------------------
def label(r):  # residue label
    return f"{r.resname}{r.resnum}"

def analyse():
    os.makedirs(OUT_DIR, exist_ok=True)
    atom_rows, res_rows = [], []
    pdbs = {}

    for ndx in sorted(glob.glob(os.path.join(DATA_DIR, "hbond_*_r*.ndx"))):
        m = re.match(r"hbond_(.+)_r(\d+)\.ndx", os.path.basename(ndx))
        system, rep = m.group(1), int(m.group(2))
        xpm = os.path.join(DATA_DIR, f"hbmap_{system}_r{rep}.xpm")
        if system not in SYSTEMS or not os.path.exists(xpm):
            print(f"skip {ndx} (system not listed or no xpm)")
            continue
        if system not in pdbs:
            pdbs[system] = read_pdb(os.path.join(DATA_DIR, f"1000ns_{system}_r1.pdb"))
        pdb = pdbs[system]

        mat = read_xpm(xpm)
        trip = read_ndx_triplets(ndx, mat.shape[0])
        occ = mat.mean(axis=1) * 100

        # atom level + residue-pair level (union of frames for same residue pair)
        pair_frames = {}
        for (d, h, a), o, row in zip(trip, occ, mat):
            D, A = pdb.loc[d], pdb.loc[a]
            dl, al = label(D), label(A)
            atom_rows.append(dict(system=system, rep=rep, donor=f"{dl}:{D.atom}",
                                  acceptor=f"{al}:{A.atom}", occupancy=o))
            d_prot, a_prot = D.resname in PROTEIN, A.resname in PROTEIN
            if d_prot != a_prot:  # protein <-> ligand/ethanol/other
                partner = al if d_prot else dl
                residue = dl if d_prot else al
                role = "protein=donor" if d_prot else "protein=acceptor"
            else:
                partner, residue, role = dl, al, "same-type"
            key = (partner, residue)
            pair_frames[key] = pair_frames.get(key, np.zeros(mat.shape[1], bool)) | row
        for (partner, residue), fr in pair_frames.items():
            res_rows.append(dict(system=system, rep=rep, partner=partner,
                                 residue=residue, occupancy=fr.mean() * 100))
        print(f"{system} r{rep}: {len(trip)} H-bonds, {mat.shape[1]} frames")

    atom_df, res_df = pd.DataFrame(atom_rows), pd.DataFrame(res_rows)
    atom_df.to_csv(f"{OUT_DIR}/hbond_atom_level.csv", index=False)
    res_df.to_csv(f"{OUT_DIR}/hbond_residue_level.csv", index=False)

    # mean +/- SD over replicates (missing = 0)
    reps = sorted(res_df.rep.unique())
    wide = res_df.pivot_table(index=["system", "partner", "residue"], columns="rep",
                              values="occupancy", fill_value=0).reindex(columns=reps, fill_value=0)
    summ = pd.DataFrame({"mean": wide.mean(axis=1), "sd": wide.std(axis=1, ddof=1),
                         "n_rep_seen": (wide > 0).sum(axis=1)}).reset_index()
    summ = pd.concat([summ, wide.reset_index(drop=True).add_prefix("r")], axis=1)
    summ["resnum"] = summ.residue.str.extract(r"(\d+)").astype(int)
    summ.sort_values(["system", "mean"], ascending=[True, False]).to_csv(
        f"{OUT_DIR}/hbond_summary.csv", index=False)
    return summ

# ----------------------------- presentation -----------------------------
def main_table(summ):
    keep = summ.groupby("residue")["mean"].max() >= MAIN_CUTOFF
    t = summ[summ.residue.isin(keep[keep].index)].copy()
    t["value"] = t["mean"].round(1).astype(str) + " +/- " + t["sd"].round(1).astype(str)
    tab = t.pivot_table(index=["partner", "residue"], columns="system", values="value",
                        aggfunc="first").reindex(columns=SYSTEMS).fillna("-")
    tab = tab.iloc[np.argsort([int(re.search(r"\d+", r[1]).group()) for r in tab.index])]
    tab.to_csv(f"{OUT_DIR}/table_main_>={int(MAIN_CUTOFF)}pct.csv")
    print(tab.to_string())

def heatmap(summ):
    keep = summ.groupby("residue")["mean"].max() >= HEATMAP_CUTOFF
    t = summ[summ.residue.isin(keep[keep].index)]
    # combine partners (e.g. ligand + ethanol) per residue by summing means
    m = t.pivot_table(index="residue", columns="system", values="mean", aggfunc="sum")
    m = m.reindex(columns=SYSTEMS).fillna(0)
    m = m.iloc[np.argsort([int(re.search(r"\d+", r).group()) for r in m.index])]
    fig, ax = plt.subplots(figsize=(1.0 * len(SYSTEMS) + 2, 0.3 * len(m) + 2))
    im = ax.imshow(m.values, aspect="auto", cmap="viridis", vmin=0, vmax=max(50, m.values.max()))
    ax.set_xticks(range(len(m.columns))); ax.set_xticklabels(m.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(m))); ax.set_yticklabels(m.index)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if v >= HEATMAP_CUTOFF:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        color="white" if v < 35 else "black", fontsize=8)
    plt.colorbar(im, label="H-bond occupancy (%)")
    ax.set_title(f"H-bond occupancy (mean of 3 replicates), residues >={int(HEATMAP_CUTOFF)}%")
    fig.tight_layout(); fig.savefig(f"{OUT_DIR}/fig_main_heatmap.png", dpi=300); plt.close(fig)

def supp_bars(summ, res_df_path):
    os.makedirs(f"{OUT_DIR}/supp", exist_ok=True)
    res_df = pd.read_csv(res_df_path)
    for s in SYSTEMS:
        d = summ[(summ.system == s) & (summ["mean"] >= SUPP_CUTOFF)]
        if d.empty:
            continue
        d = d.assign(lab=d.partner + "-" + d.residue).sort_values("resnum")
        fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(d) + 2), 4))
        ax.bar(d.lab, d["mean"], yerr=d.sd, capsize=3, color="steelblue", alpha=.8)
        for i, (_, r) in enumerate(d.iterrows()):   # individual replicates as dots
            pts = res_df[(res_df.system == s) & (res_df.partner == r.partner) &
                         (res_df.residue == r.residue)].occupancy
            ax.scatter([i] * len(pts), pts, color="k", s=10, zorder=3)
        ax.set_ylabel("H-bond occupancy (%)"); ax.set_title(f"{s} (mean +/- SD, dots = replicates)")
        plt.xticks(rotation=60, ha="right"); fig.tight_layout()
        fig.savefig(f"{OUT_DIR}/supp/bar_{s}.png", dpi=300); plt.close(fig)

if __name__ == "__main__":
    summ = analyse()
    main_table(summ)
    heatmap(summ)
    supp_bars(summ, f"{OUT_DIR}/hbond_residue_level.csv")
    print(f"\nDone. See {OUT_DIR}/")