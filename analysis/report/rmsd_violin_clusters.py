"""RMSD to the crystal ligand per barrier-cut cluster of each complex (Shu-Yu): one panel per complex,
one violin per cluster of that complex's own barrier-cut clustering (report_cut.sbatch, kT_cut 1,
<work>/<c>/<label>-cut/labels.csv), clusters >= --min-mass sorted by population, plus "all frames".
RMSD as rmsd_violin.py / rmsd_density.py (heavy atoms superposed on the crystal ligand, best over the
symmetry-equivalent graph matches). The gold cluster is the one the crystal falls in.

Also writes, for a look in PyMOL, the solute (protein + ligand, no water) of the --n-worst frames with
the largest RMSD and of each cluster's central frame (cluster_cut's groups.csv):
  <work>/<c>/rmsd_clusters/<label>_worst_<rank>_frame<i>.pdb, <label>_central_<cluster>.pdb

  python -I rmsd_violin_clusters.py --name folic-acid  -> results/<name>/<version>/<name>_rmsd_violin_clusters.png/.csv
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rmsd_density import rmsd_to_crystal  # noqa: E402

WORK = "/scicore/home/lill0000/popa0002/thesis/analysis-report"
REPO = os.path.dirname(os.path.dirname(HERE))
GOLD, COL, ALL = "#e6a100", "#b2182b", "#9a9a9a"


def solute_frames(c, frames, path_for):
    """Write the full solute (protein + ligand) of the given complex frames."""
    run = c["run"]
    t = md.load(f"{run}/REST2-run1/solute_state0_prod1.nc", top=f"{run}/build/built.solute.pdb")
    for i, p in zip(frames, path_for):
        t[int(i)].save_pdb(p)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True)
    ap.add_argument("--min-mass", type=float, default=0.01)
    ap.add_argument("--n-worst", type=int, default=3)
    args = ap.parse_args()
    man = json.load(open(f"{WORK}/{args.name}/manifest.json"))
    out_dir = f"{WORK}/{args.name}/rmsd_clusters"
    os.makedirs(out_dir, exist_ok=True)
    panels, rows = [], []
    for c in man["complexes"]:
        label = c["label"]
        rmsd = rmsd_to_crystal(c, c["refs"][0])
        cut = f"{WORK}/{args.name}/{label}-cut"
        lab = pd.read_csv(f"{cut}/labels.csv", keep_default_na=False, dtype={"cluster": str})["cluster"].to_numpy()
        if len(lab) != len(rmsd):
            sys.exit(f"{label}: {len(lab)} labels, {len(rmsd)} frames")
        crystal = {v["cluster"] for v in json.load(open(f"{cut}/summary.json")).get("crystal_assignment", {}).values()
                   if v.get("cluster")}
        names, counts = np.unique(lab, return_counts=True)
        shown = [(n, k / len(lab)) for n, k in sorted(zip(names, counts), key=lambda t: -t[1]) if k / len(lab) >= args.min_mass]
        rest = 1 - sum(m for _, m in shown)
        panels.append((label, rmsd, lab, shown, rest, crystal))
        groups = pd.read_csv(f"{cut}/groups.csv", keep_default_na=False, dtype={"cluster": str}).set_index("cluster")
        for n, m in shown:
            v = rmsd[lab == n]
            rows.append({"complex": label, "cluster": n, "mass": round(m, 4), "median_rmsd": round(float(np.median(v)), 3),
                         "q25": round(float(np.percentile(v, 25)), 3), "q75": round(float(np.percentile(v, 75)), 3),
                         "crystal_cluster": n in crystal, "central_frame": int(groups.loc[n, "central_frame"])})
        worst = np.argsort(rmsd)[::-1][:args.n_worst]
        solute_frames(c, list(worst) + [int(groups.loc[n, "central_frame"]) for n, _ in shown],
                      [f"{out_dir}/{label}_worst_{r + 1}_frame{i}.pdb" for r, i in enumerate(worst)]
                      + [f"{out_dir}/{label}_central_{n}.pdb" for n, _ in shown])
        print(f"{args.name} {label}: worst frames " + ", ".join(f"{i} ({rmsd[i]:.2f} A, cluster {lab[i]})" for i in worst)
              + " | per cluster median: " + ", ".join(f"{r['cluster']} {r['median_rmsd']:.2f}" for r in rows if r["complex"] == label),
              flush=True)
    stem = os.path.join(REPO, "results", args.name, man["version"], f"{args.name}_rmsd_violin_clusters")
    pd.DataFrame(rows).to_csv(stem + ".csv", index=False)

    widths = [1 + len(p[3]) for p in panels]
    fig, axes = plt.subplots(1, len(panels), figsize=(1.5 + 0.75 * sum(widths), 5.4), layout="constrained",
                             gridspec_kw={"width_ratios": widths}, squeeze=False, sharey=True)
    for ax, (label, rmsd, lab, shown, rest, crystal) in zip(axes[0], panels):
        series = [("all", 1.0, rmsd, ALL)] + [(n, m, rmsd[lab == n], GOLD if n in crystal else COL) for n, m in shown]
        for x, (n, m, v, col) in enumerate(series):
            vp = ax.violinplot(v, positions=[x], widths=0.8, showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor(col); b.set_edgecolor("k"); b.set_linewidth(0.7); b.set_alpha(0.85)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.vlines(x, q1, q3, color="k", lw=3.5, zorder=3)
            ax.scatter([x], [med], color="white", edgecolor="k", s=22, zorder=4)
        ax.set_xticks(range(len(series)), ["all\nframes"] + [f"{n}\n{100 * m:.1f}%" for n, _, _, _ in series[1:]
                                                            for m in [dict(shown)[n]]], fontsize=8.5)
        for t, (n, *_ ) in zip(ax.get_xticklabels()[1:], series[1:]):
            if n in crystal:
                t.set_color(GOLD); t.set_fontweight("bold")
        ax.set_title(f"{label} complex: {len(shown)} clusters ≥ {100 * args.min_mass:g} %"
                     + (f"\n(the rest hold {100 * rest:.1f} %)" if rest > 0.0005 else ""), fontsize=10)
        ax.set_xlim(-0.6, len(series) - 0.4)
        ax.grid(alpha=0.3, axis="y")
    axes[0][0].set_ylabel("heavy-atom RMSD to that complex's crystal ligand (Å)")
    axes[0][0].set_ylim(bottom=0)
    fig.suptitle(f"{args.name}: ligand conformation vs the crystal, per cluster of the complex (barrier cuts, kT_cut 1)\n"
                 "REST2 state 0; ligand superposed on the crystal ligand, symmetry-corrected; dot = median, bar = IQR; "
                 "gold = cluster the crystal falls in", fontsize=10)
    fig.savefig(stem + ".png", dpi=160)
    print(f"wrote {stem}.png / .csv; structures in {out_dir}")


if __name__ == "__main__":
    main()
