"""Conformational entropy vs the barrier-cut resolution: S_c = -sum pi_c ln pi_c (k_B) of the torsion-
partitioning clusters (cut_clustering, MD-tools 306b109) over a grid of the two settings,
  pi_min = mass_floor (smallest piece, fraction of all frames)   and   E* = kT_cut (lowest barrier cut),
for the water run and every complex of one compound (two panels: water, complexes). Same inputs as
report_cut.sbatch (prepare.py's trajectories and torsions, every frame); the molecule's symmetry is
enumerated once per system and reused for every setting. Shaded: pi_min below 10 frames of that run
(clusters that small are sampling noise). Error bars: block SE over --blocks time blocks of the labels.

  python -I entropy_scan.py --name folic-acid   -> results/<name>/<version>/<name>_entropy_scan.png/.csv
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "shuyu_clustering"))
from kl_cross import WORK, REPO, KINDS, thdbscan  # noqa: E402
from cluster_prune import load_torsions  # noqa: E402

FLOORS = (1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4)
KTS = (1.0, 2.0, 3.0)
COLOURS = {1.0: "#1b9e77", 2.0: "#d95f02", 3.0: "#7570b3"}


def entropy(labels):
    _, n = np.unique(labels, return_counts=True)
    p = n / n.sum()
    return float(-(p * np.log(p)).sum())


def scan(system, extra, blocks):
    from md_tools.analysis import cut_clustering
    tors = load_torsions([f"{system['torsions']}/{k}.txt" for k in KINDS
                          if os.path.getsize(f"{system['torsions']}/{k}.txt") > 0])
    t = md.load(system["traj"], top=system["top"])
    mol = Chem.SDMolSupplier(system["sdf"], removeHs=False)[0]
    fit = thdbscan(t, mol, list(tors), list(tors.values()), extra)
    rows = []
    for kt in KTS:
        for floor in FLOORS:
            r = cut_clustering(fit, kT_cut=kt, mass_floor=floor)
            lab = np.asarray(r["labels"])
            se = np.std([entropy(b) for b in np.array_split(lab, blocks)], ddof=1) / np.sqrt(blocks)
            rows.append({"E_star": kt, "pi_min": floor, "n_clusters": int(r["n_clusters"]),
                         "S": entropy(lab), "S_se": float(se), "frames": len(lab)})
            print(f"  {system['label']:6s} E* {kt:g}  pi_min {floor:g}: {r['n_clusters']:4d} clusters, "
                  f"S {rows[-1]['S']:.3f} +- {se:.3f}", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True)
    ap.add_argument("--blocks", type=int, default=5)
    args = ap.parse_args()
    man = json.load(open(f"{WORK}/{args.name}/manifest.json"))
    extra = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
    res = {}
    for s in [man["water"]] + man["complexes"]:
        res[s["label"]] = scan(s, extra, args.blocks)
    df = pd.concat([d.assign(system=k) for k, d in res.items()], ignore_index=True)
    out = os.path.join(REPO, "results", args.name, man["version"], f"{args.name}_entropy_scan")
    df.to_csv(out + ".csv", index=False)

    from matplotlib.lines import Line2D
    cplx = [c["label"] for c in man["complexes"]]
    styles = [("-", "o"), ("--", "s"), (":", "^"), ("-.", "D")]       # one line style + marker per complex
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), layout="constrained")
    for ax, (title, systems) in zip(axes, [("water", ["water"]), ("complex", cplx)]):
        for k, s in enumerate(systems):
            ls, mk = styles[k % len(styles)] if s != "water" else ("-", "o")
            for kt in KTS:
                d = res[s][res[s].E_star == kt].sort_values("pi_min")
                ax.errorbar(d.pi_min, d.S, yerr=d.S_se, color=COLOURS[kt], ls=ls, marker=mk, ms=5, capsize=2, lw=1.6)
                if kt == 1.0:                       # the default setting used everywhere else
                    i = int(np.argmin(np.abs(np.log(d.pi_min.to_numpy()) - np.log(1e-3))))
                    ax.plot(d.pi_min.iloc[i], d.S.iloc[i], "*", color="k", ms=14, zorder=5)
            n = int(res[s].frames.iloc[0])
            ax.axvspan(1e-4 * 0.8, 10.0 / n, color="0.88", zorder=0, lw=0)
        handles = [Line2D([], [], color=COLOURS[kt], lw=2.5, label=f"E* = {kt:g} k$_B$T") for kt in KTS]
        if len(systems) > 1:
            handles += [Line2D([], [], color="k", ls=styles[k][0], marker=styles[k][1], ms=5, lw=1.6,
                               label=f"{s} complex") for k, s in enumerate(systems)]
        handles.append(Line2D([], [], color="k", marker="*", ls="", ms=12, label="default (E* 1, π$_{min}$ 10$^{-3}$)"))
        ax.legend(handles=handles, fontsize=8.5, handlelength=3.2, loc="lower right" if len(systems) == 1 else "upper right")
        ax.set_xscale("log")
        ax.set_xlim(1.25e-2, 0.8e-4)                # coarse (left) -> fine (right)
        ax.set_xlabel("π$_{min}$ (mass floor, fraction of all frames)")
        ax.set_title(title if len(systems) == 1 else f"complex ({' / '.join(systems)})", fontsize=11)
        ax.grid(alpha=0.3, which="both")
    axes[0].set_ylabel("S$_c$ = −Σ π$_c$ ln π$_c$  (k$_B$)")
    top = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, top)
    fig.suptitle(f"{args.name}: conformational entropy vs barrier-cut resolution (torsion partitioning, REST2 state 0, every frame)\n"
                 "π$_c$ = fraction of the run's frames in cluster c; grey = π$_{min}$ below 10 frames of that run; "
                 f"bars: block SE over {args.blocks} time blocks", fontsize=11)
    fig.savefig(out + ".png", dpi=150)
    print(f"wrote {out}.png / .csv")


if __name__ == "__main__":
    main()
