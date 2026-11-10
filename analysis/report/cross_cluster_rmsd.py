"""Cluster-by-cluster RMSD between the bound and the free ensemble (Shu-Yu): rows = the complex's own
barrier-cut clusters, columns = the water run's own barrier-cut clusters (two separate clusterings, the
report_cut.sbatch runs <work>/<c>/<system>-cut), each axis sorted by population, largest first.

  cell (X^p, Y^w) = mean symmetry-corrected heavy-atom RMSD (superposed, minimised over the molecule's
                    symmetry operations) between complex frames of X and water frames of Y (A)
A low cell: that bound conformer family is, in 3D, one the free ligand also visits; a row with no low cell:
the bound conformers are far from everything water samples. Two clusterings, so the matrix is not
symmetric and its diagonal is not special. Clusters below --min-mass (default 1 %) of their own run are
left out (titles say how much); at most --max-per-cluster random frames per cluster (seeded). The cluster
each crystal falls in is marked gold on both axes (cluster_cut's crystal assignment).
Both ensembles in the water ligand's atom order (kl_cross.ensembles), so the RMSD compares the same atoms.

  python -I cross_cluster_rmsd.py --name folic-acid     -> results/<name>/<version>/<name>_<label>_cut_cross_rmsd.png/.csv
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kl_cross import WORK, REPO, ensembles  # noqa: E402
from md_tools.analysis import enumerate_symmetry, heavy_atoms, operations_on_atoms, symmetric_rmsd  # noqa: E402

GOLD = "#e6a100"


def clusters(cut_dir, n_frames, min_mass):
    """(labels per frame, [(name, mass)] >= min_mass sorted by size, left-out mass, crystal clusters)."""
    s = json.load(open(f"{cut_dir}/summary.json"))
    lab = pd.read_csv(f"{cut_dir}/labels.csv", keep_default_na=False, dtype={"cluster": str})["cluster"].to_numpy()
    if int(s["stride"]) != 1 or len(lab) != n_frames:
        sys.exit(f"{cut_dir}: {len(lab)} labels (stride {s['stride']}) for {n_frames} frames")
    names, cnt = np.unique(lab, return_counts=True)
    order = [(str(n), c / len(lab)) for n, c in sorted(zip(names, cnt), key=lambda t: -t[1])]
    shown = [(n, m) for n, m in order if m >= min_mass]
    crystal = {v["cluster"] for v in s.get("crystal_assignment", {}).values() if v.get("cluster")}
    return lab, shown, 1.0 - sum(m for _, m in shown), len(order), crystal


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True)
    ap.add_argument("--min-mass", type=float, default=0.01)
    ap.add_argument("--max-per-cluster", type=int, default=400)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()
    man = json.load(open(f"{WORK}/{args.name}/manifest.json"))
    mol = Chem.SDMolSupplier(man["water"]["sdf"], removeHs=False)[0]
    heavy = heavy_atoms(mol)
    kw = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
    ops = operations_on_atoms(enumerate_symmetry(mol, **kw), heavy)
    rng = np.random.default_rng(0)
    outdir = os.path.join(REPO, "results", args.name, man["version"])
    for c in man["complexes"]:
        tw, tc, _, _ = ensembles(man, c)
        Xw = tw.xyz[:, heavy, :] * 10.0                   # both in the water ligand's atom order, A
        Xp = tc.xyz[:, heavy, :] * 10.0
        lw, cw, rest_w, nw, xw = clusters(f"{WORK}/{args.name}/water-cut", tw.n_frames, args.min_mass)
        lp, cp, rest_p, np_, xp = clusters(f"{WORK}/{args.name}/{c['label']}-cut", tc.n_frames, args.min_mass)
        pick = lambda lab, n: (lambda idx: idx if len(idx) <= args.max_per_cluster else
                               np.sort(rng.choice(idx, args.max_per_cluster, replace=False)))(np.flatnonzero(lab == n))
        W = {n: Xw[pick(lw, n)] for n, _ in cw}
        P = {n: Xp[pick(lp, n)] for n, _ in cp}
        mean = np.zeros((len(cp), len(cw)))
        std = np.zeros_like(mean)
        for i, (a, _) in enumerate(cp):
            for j, (b, _) in enumerate(cw):
                d = symmetric_rmsd(P[a], W[b], ops, workers=args.workers)
                mean[i, j], std[i, j] = d.mean(), d.std()
        stem = f"{outdir}/{args.name}_{c['label']}_cut_cross_rmsd"
        pd.DataFrame([{"complex_cluster": a, "complex_mass": round(ma, 5), "water_cluster": b, "water_mass": round(mb, 5),
                       "rmsd_mean": round(mean[i, j], 4), "rmsd_std": round(std[i, j], 4)}
                      for i, (a, ma) in enumerate(cp) for j, (b, mb) in enumerate(cw)]).to_csv(stem + ".csv", index=False)

        nr, nc = mean.shape
        fig, ax = plt.subplots(figsize=(3.5 + 0.62 * nc, 2.6 + 0.5 * nr), layout="constrained")
        im = ax.imshow(mean, cmap="viridis", vmin=0, aspect="auto")
        if nc <= 25:
            for i in range(nr):
                for j in range(nc):
                    ax.text(j, i, f"{mean[i, j]:.1f}", ha="center", va="center", fontsize=7.5 if nc > 12 else 9,
                            color="w" if mean[i, j] < 0.6 * mean.max() else "k")
        ax.set_xticks(range(nc), [f"{b}$^w$\n{100 * m:.1f}%" for b, m in cw], fontsize=8)
        ax.set_yticks(range(nr), [f"{a}$^p$  {100 * m:.1f}%" for a, m in cp], fontsize=8.5)
        for t, (b, _) in zip(ax.get_xticklabels(), cw):
            if b in xw:
                t.set_color(GOLD); t.set_fontweight("bold")
        for t, (a, _) in zip(ax.get_yticklabels(), cp):
            if a in xp:
                t.set_color(GOLD); t.set_fontweight("bold")
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_xlabel(f"water clusters (barrier cuts), largest first: {nc} of {nw} shown, the rest hold {100 * rest_w:.1f} %")
        ax.set_ylabel(f"{c['label']} complex clusters, largest first\n{nr} of {np_} shown, the rest hold {100 * rest_p:.1f} %")
        fig.colorbar(im, ax=ax, shrink=0.8, label="mean heavy-atom RMSD (Å), symmetry-corrected")
        fig.suptitle(f"{args.name}: bound vs free conformers, cluster by cluster ({c['label']} complex frames x water frames)\n"
                     f"gold = cluster the crystal falls in; at most {args.max_per_cluster} frames per cluster; "
                     f"clusters >= {100 * args.min_mass:g} % of their run", fontsize=10)
        fig.savefig(stem + ".png", dpi=150)
        plt.close(fig)
        best = mean.min(axis=1)
        print(f"{args.name} {c['label']}: {nr} complex x {nc} water clusters; closest water cluster per complex cluster "
              + ", ".join(f"{a} {best[i]:.2f}" for i, (a, _) in enumerate(cp)) + f" A -> {stem}.png", flush=True)


if __name__ == "__main__":
    main()
