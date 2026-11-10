"""Ligand conformation RMSD to each crystal structure, water vs complex, as violins: one figure per compound.

x: the compound's crystal structures (one per complex, e.g. folate 4LRH and 3QL3); at each, two violins:
the ligand in water and in that complex, RMSD of every frame to that crystal ligand. RMSD as in
rmsd_density.py: heavy atoms superposed on the crystal (Kabsch), lowest over the symmetry-equivalent graph
matches -- the ligand's SHAPE, not its place in the pocket. REST2 state 0, every saved frame. Inside each
violin: median (white dot) and interquartile range (bar).

  python -I rmsd_violin.py --name folic-acid     -> results/<name>/<version>/<name>_rmsd_violin.png
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rmsd_density import rmsd_to_crystal  # noqa: E402

WORK = "/scicore/home/lill0000/popa0002/thesis/analysis-report"
REPO = os.path.dirname(os.path.dirname(HERE))
WATER, COMPLEX = "#9ecae1", "#b2182b"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True)
    args = ap.parse_args()
    man = json.load(open(f"{WORK}/{args.name}/manifest.json"))
    data = []                                    # (crystal label, water rmsd, complex rmsd)
    for c in man["complexes"]:
        ref = c["refs"][0]
        data.append((c["label"], rmsd_to_crystal(man["water"], ref), rmsd_to_crystal(c, ref)))
        print(f"{args.name} {c['label']}: median RMSD water {np.median(data[-1][1]):.2f} A, "
              f"complex {np.median(data[-1][2]):.2f} A", flush=True)

    fig, ax = plt.subplots(figsize=(2.2 + 2.4 * len(data), 5.2), layout="constrained")
    w = 0.38
    for i, (label, rw, rc) in enumerate(data):
        for x, vals, col in ((i - w / 2 - 0.02, rw, WATER), (i + w / 2 + 0.02, rc, COMPLEX)):
            v = ax.violinplot(vals, positions=[x], widths=w, showextrema=False)
            for b in v["bodies"]:
                b.set_facecolor(col); b.set_edgecolor("k"); b.set_linewidth(0.8); b.set_alpha(0.85)
            q1, med, q3 = np.percentile(vals, [25, 50, 75])
            ax.vlines(x, q1, q3, color="k", lw=4, zorder=3)
            ax.scatter([x], [med], color="white", edgecolor="k", s=28, zorder=4)
    ax.set_xticks(range(len(data)), [f"crystal {d[0]}" for d in data])
    ax.set_xlim(-0.6, len(data) - 0.4)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("heavy-atom RMSD to the crystal ligand (Å)")
    ax.grid(alpha=0.3, axis="y")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=WATER, edgecolor="k", label="in water"),
                       Patch(facecolor=COMPLEX, edgecolor="k", label="in the complex")],
              loc="upper left", fontsize=9)
    ax.set_title(f"{args.name}: ligand conformation vs the crystal structure\nREST2 state 0, every frame; "
                 "ligand superposed on the crystal ligand,\nsymmetry-corrected; dot = median, bar = interquartile range",
                 fontsize=10)
    out = os.path.join(REPO, "results", args.name, man["version"], f"{args.name}_rmsd_violin.png")
    fig.savefig(out, dpi=160)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
