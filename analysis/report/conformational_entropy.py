"""Conformational entropy from cluster populations, as a flexibility measure comparable across ligands.

  S_c = - sum_c pi_c ln pi_c   (units of k_B),   e^S_c = the effective number of equally populated
                                                conformers ("behaves as if it had 1.84 conformers")

for every ligand in water and in each complex, from the clustering outputs written by report.sbatch
(symmetry-first HDBSCAN, --no-single-cluster, MD-tools c30c1c0) and report_cut.sbatch (torsion
partitioning by barrier cuts, MD-tools 306b109), one figure per clustering method. Points in one
figure come from exactly one algorithm with one set of settings: the script reads each run's own
record (summary.json) and stops if route, settings or md-tools version differ between points.

  pi_c   HDBSCAN: cluster masses renormalised over the frames assigned to a cluster (noise and
         ambiguous frames are left out; their fraction is in the CSV). A run where no cluster is found
         (NO_CLUSTER.txt: one basin, e.g. 4ibi bound) is drawn at e^S = 1 with a hollow marker.
         Barrier cuts: every frame is labelled; all clusters.
  error  block standard error of e^S over --blocks consecutive time blocks of the state-0 labels
  x      rotatable bonds, RDKit on the heavy-atom graph (the definition used to select the ligands)

  python -I conformational_entropy.py [--compounds folic-acid paracetamol 4ibi ...] [--out DIR]
writes conformational_entropy_hdbscan.png, conformational_entropy_cut.png and conformational_entropy.csv.
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

WORK = "/scicore/home/lill0000/popa0002/thesis/analysis-report"
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import csv as _csv
DEFAULT = [r["name"] for r in _csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "compounds.csv")))]
METHODS = {"hdbscan": ("clustering", "symmetry-first HDBSCAN (no single cluster)"),
           "cut": ("cut", "torsion partitioning (barrier cuts)")}


def entropy(p):
    p = np.asarray([x for x in p if x > 0], float)
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def letters(i):
    """0 -> A, 25 -> Z, 26 -> AA: MD-tools' group_name, the names labels.csv uses."""
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def block_se(labels, clusters, blocks):
    """Standard error of e^S over consecutive blocks of the per-frame labels (unassigned dropped)."""
    vals = []
    for b in np.array_split(np.asarray(labels), blocks):
        b = b[np.isin(b, clusters)]
        if len(b):
            _, counts = np.unique(b, return_counts=True)
            vals.append(np.exp(entropy(counts)))
    return float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else float("nan")


def read_run(d, method, blocks):
    """(populations, unassigned fraction, e^S block SE, algorithm record, note) of one clustering output."""
    if method == "hdbscan" and os.path.exists(f"{d}/NO_CLUSTER.txt"):
        return [1.0], 1.0, 0.0, None, "no cluster found (one basin): drawn at e^S = 1"
    s = json.load(open(f"{d}/summary.json"))
    lab = pd.read_csv(f"{d}/labels.csv").iloc[:, 2].astype(str).to_numpy()
    if method == "hdbscan":
        masses = s["cluster_mass"]
        p = list(masses.values()) if isinstance(masses, dict) else list(masses)
        names = [letters(i) for i in range(len(p))]          # labels.csv names clusters A, B, ... by index
        rec = {"route": s["route"], "allow_single_cluster": s["allow_single_cluster"], "mass_floor": s["mass_floor"],
               "vote": json.dumps(s.get("vote"), sort_keys=True, default=str)}
        unassigned = float(s.get("unassigned_mass", 1.0 - sum(p)))
        clusters = [str(n) for n in names]
    else:
        p = list(s["cluster_mass"].values())
        st = {k: v for k, v in s["settings"].items() if k not in ("floor_frames", "reference")}
        rec = {"method": s["method"], "md_tools": s["md_tools"], **{k: json.dumps(v) for k, v in st.items()}}
        unassigned = 0.0
        clusters = list(s["cluster_mass"])
    return p, unassigned, block_se(lab, clusters, blocks), rec, ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--compounds", nargs="+", default=DEFAULT)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", default=os.path.join(REPO, "results", "overview", "conformational-entropy"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rows = []
    for name in args.compounds:
        man = json.load(open(f"{WORK}/{name}/manifest.json"))
        heavy = Chem.RemoveHs(Chem.MolFromMolFile(man["water"]["sdf"]))
        rb = rdMolDescriptors.CalcNumRotatableBonds(heavy)
        systems = [("water", "water")] + [(c["label"], f"{c['label']} complex") for c in man["complexes"]]
        for label, system in systems:
            for method, (suffix, _) in METHODS.items():
                d = f"{WORK}/{name}/{label}-{suffix}"
                if not os.path.isdir(d):
                    sys.exit(f"{d}: missing -- run report.sbatch / report_cut.sbatch for {name}")
                p, unassigned, se, rec, note = read_run(d, method, args.blocks)
                S = entropy(p)
                rows.append({"method": method, "ligand": name, "system": system, "rot_bonds": rb,
                             "n_clusters": len(p) if not note else 0, "S_kB": round(S, 4),
                             "exp_S": round(float(np.exp(S)), 3), "exp_S_block_se": round(se, 3),
                             "unassigned_fraction": round(unassigned, 4), "note": note,
                             "algorithm": json.dumps(rec, sort_keys=True) if rec else ""})
    df = pd.DataFrame(rows)

    # one algorithm, one set of settings per figure
    for method in METHODS:
        recs = set(df[(df.method == method) & (df.algorithm != "")].algorithm)
        if len(recs) != 1:
            sys.exit(f"{method}: the clustering outputs were made with different settings:\n  " + "\n  ".join(recs))
        print(f"{method}: one algorithm record for all points: {recs.pop()}")
    df.drop(columns="algorithm").to_csv(f"{args.out}/conformational_entropy.csv", index=False)
    with pd.option_context("display.width", 200):
        print(df.drop(columns="algorithm").to_string(index=False))

    palette = list(plt.get_cmap("tab20").colors[::2]) + list(plt.get_cmap("tab20").colors[1::2]) + list(plt.get_cmap("Dark2").colors)
    colours = dict(zip(args.compounds, palette))
    markers = ["s", "D", "^", "v"]
    for method, (_, title) in METHODS.items():
        sub = df[df.method == method]
        fig, ax = plt.subplots(figsize=(13.5, 6.8), layout="constrained")
        per_rb = {}
        for name in args.compounds:                       # spread ligands that share a bond count
            rb = int(sub[sub.ligand == name].rot_bonds.iloc[0])
            per_rb.setdefault(rb, []).append(name)
        offset = {n: (i - (len(v) - 1) / 2) * min(0.18, 0.8 / max(len(v), 1)) for v in per_rb.values() for i, n in enumerate(v)}
        for name in args.compounds:
            g = sub[sub.ligand == name]
            col = colours[name]
            x0 = g.rot_bonds.iloc[0] + offset[name]
            w = g[g.system == "water"].iloc[0]
            cpx = g[g.system != "water"]
            for k, (_, c) in enumerate(cpx.iterrows()):
                ax.plot([x0 - 0.05, x0 + 0.05], [w.exp_S, c.exp_S], color=col, lw=0.8, ls=":", alpha=0.7)
            ax.errorbar(x0 - 0.05, w.exp_S, yerr=w.exp_S_block_se, fmt="o", mfc="white", mec=col, ecolor=col,
                        ms=9, mew=2, capsize=3, label=f"{name} ({int(w.rot_bonds)} rot. bonds)")
            for k, (_, c) in enumerate(cpx.iterrows()):
                hollow = bool(c.note)
                ax.errorbar(x0 + 0.05, c.exp_S, yerr=None if hollow else c.exp_S_block_se,
                            fmt=markers[k % len(markers)], mfc="white" if hollow else col, mec=col, ecolor=col,
                            ms=9, mew=2, capsize=3)
                ax.annotate(c.system.replace(" complex", ""), (x0 + 0.05, c.exp_S), textcoords="offset points",
                            xytext=(7, -3), fontsize=7.5, color=col)
        ax.set_yscale("log")
        ax.set_xlabel("rotatable bonds (heavy-atom graph)")
        ax.set_ylabel(r"$e^{S_c}$  (effective number of equally populated conformers)")
        ax.set_xticks(sorted(df.rot_bonds.unique()))
        ax.grid(alpha=0.3, which="both")
        h = [plt.Line2D([], [], marker="o", ls="", mfc="white", mec="0.3", mew=2, ms=8, label="in water"),
             plt.Line2D([], [], marker="s", ls="", color="0.3", ms=8, label="in its complex (PDB id at the point)")]
        if (sub.note != "").any():
            h.append(plt.Line2D([], [], marker="s", ls="", mfc="white", mec="0.3", mew=2, ms=8,
                                label="complex, no cluster found: one basin, e^S = 1"))
        leg1 = ax.legend(handles=h, loc="upper left", fontsize=8.5)
        ax.add_artist(leg1)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, title="ligand",
                  ncol=1 if len(args.compounds) <= 12 else 2)
        extra = ("populations over the assigned frames (noise and ambiguous frames left out)"
                 if method == "hdbscan" else "every frame labelled")
        ax.set_title(f"Conformational entropy $S_c = -\\sum_c \\pi_c \\ln \\pi_c$ from cluster populations\n"
                     f"{title}, REST2 state 0\n{extra}; bars: block SE over {args.blocks} time blocks", fontsize=10)
        fig.savefig(f"{args.out}/conformational_entropy_{method}.png", dpi=150)
        plt.close(fig)
    print(f"wrote {args.out}/conformational_entropy_hdbscan.png, _cut.png and conformational_entropy.csv")


if __name__ == "__main__":
    main()
