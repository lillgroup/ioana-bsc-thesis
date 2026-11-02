"""Torsions of the ligand in water vs in its complex(es): one figure per compound, one panel per
torsion, the state-0 distributions overlaid (all frames, no clustering) and each complex's crystal
value as a dashed line in that complex's colour.

Generalises ../shuyu_clustering/compare_water_complex.py (folate only) to any compound in the
manifest (prepare.py). Each system reads its own torsion files, which prepare.py mapped from the water
ligand, so t<k> is the same bond here, in the clustering figures and in every other plot.

Reading a panel: a complex peak on a water peak or tail -> the protein SELECTS a conformer water
already samples; a complex peak where water has ~nothing -> the protein INDUCES it; a complex peak
away from its own dashed line -> the complex does not hold the crystal value.

  python -I water_vs_complex.py --manifest <work>/<compound>/manifest.json --out <prefix>
writes <prefix>.png and <prefix>.csv (per torsion x complex: peak, water mass at that peak,
complex and water mass at the crystal value, +-30 deg).
"""
import argparse
import io
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "shuyu_clustering"))
from cluster_prune import load_torsions, reference_torsions, molecule_png, simulated_ns  # noqa: E402

KINDS = ("torsions", "torsions_rigid", "torsions_polarH")
COLOURS = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e"]


def wrap(x):
    return np.angle(np.exp(1j * np.asarray(x)))


def mass_within(values, centre, window):
    return float(np.mean(np.abs(wrap(values - centre)) <= window))


def peak(values, edges):
    h, _ = np.histogram(values, bins=edges)
    return 0.5 * (edges[h.argmax()] + edges[h.argmax() + 1])


def system_torsions(system, stride):
    """(names, per-frame dihedrals (F, n), quads, RDKit mol, simulated ns) for one system."""
    tors = load_torsions([f"{system['torsions']}/{k}.txt" for k in KINDS
                          if os.path.getsize(f"{system['torsions']}/{k}.txt") > 0])
    quads = list(tors.values())
    t = md.load(system["traj"], top=system["top"], stride=stride)
    mol = Chem.SDMolSupplier(system["sdf"], removeHs=False)[0]
    return list(tors), md.compute_dihedrals(t, quads), quads, mol, simulated_ns(t, stride)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--bins", type=int, default=72, help="histogram bins over 360 deg")
    ap.add_argument("--window", type=float, default=30.0, help="+-deg for the populations in the CSV")
    args = ap.parse_args()
    man = json.load(open(args.manifest))
    stride = man["stride"]

    names, wvals, wquads, wmol, wns = system_torsions(man["water"], stride)
    systems = {"water": (wvals, "#9a9a9a", None, wns)}
    for k, c in enumerate(man["complexes"]):
        cnames, cvals, cquads, cmol, cns = system_torsions(c, stride)
        assert len(cnames) == len(names), f"{c['label']}: {len(cnames)} torsions, water has {len(names)}"
        ref = reference_torsions(c["refs"][0], cmol, cquads)
        systems[f"{c['label']} complex"] = (cvals, COLOURS[k % len(COLOURS)], ref, cns)

    edges = np.linspace(-np.pi, np.pi, args.bins + 1)
    win = np.deg2rad(args.window)
    atom_names = [[wmol.GetAtomWithIdx(i).GetSymbol() + str(i) for i in q] for q in wquads]
    rows = []
    for k in range(len(names)):
        for label, (vals, _, ref, _) in systems.items():
            if ref is None:
                continue
            p = peak(vals[:, k], edges)
            row = {"torsion": f"t{k}", "atoms": "-".join(atom_names[k]), "complex": label,
                   "complex_peak_deg": round(float(np.rad2deg(p)), 1),
                   "complex_mass_at_peak": round(mass_within(vals[:, k], p, win), 3),
                   "water_mass_at_complex_peak": round(mass_within(wvals[:, k], p, win), 3),
                   "crystal_deg": ", ".join(f"{np.rad2deg(v):.1f}" for v in ref[k]) or "n/a"}
            if ref[k]:
                row["complex_mass_at_crystal"] = round(max(mass_within(vals[:, k], v, win) for v in ref[k]), 3)
                row["water_mass_at_crystal"] = round(max(mass_within(wvals[:, k], v, win) for v in ref[k]), 3)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(args.out + ".csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(df.to_string(index=False))

    ntor = len(names)
    ncol = min(4, ntor)
    nrow = -(-ntor // ncol)
    panel_w, panel_h, mol_w = 2.7, 2.1, 5.0
    fig = plt.figure(figsize=(mol_w + ncol * panel_w, max(nrow * panel_h, 4.0) + 1.2), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=[mol_w, ncol * panel_w])
    mol_ax = fig.add_subplot(outer[0])
    png = molecule_png(wmol, wquads, ["#333333"] * ntor, size=(1100, 1100))
    mol_ax.imshow(mpimg.imread(io.BytesIO(png), format="png"))
    mol_ax.axis("off")
    mol_ax.set_title("t<k> marks the central bond of torsion k", fontsize=10)
    grid = outer[1].subgridspec(nrow, ncol)
    xt = np.deg2rad([-180, -90, 0, 90, 180])
    for k, name in enumerate(names):
        ax = fig.add_subplot(grid[divmod(k, ncol)])
        for label, (vals, colour, ref, _) in systems.items():
            if label == "water":
                ax.hist(vals[:, k], bins=edges, density=True, histtype="stepfilled", color=colour, alpha=0.55,
                        lw=0, label=label)
            else:
                ax.hist(vals[:, k], bins=edges, density=True, histtype="step", color=colour, lw=1.6, label=label)
                for n, v in enumerate(ref[k]):
                    ax.axvline(v, color=colour, ls="--", lw=1.3, label=(f"crystal {label.split()[0]}" if n == 0 else None))
        ax.set_xlim(-np.pi, np.pi)
        ax.set_xticks(xt, ["-180", "-90", "0", "90", "180"])
        ax.set_yticks([])
        ax.tick_params(labelsize=8)
        kind = name.rsplit("_", 1)[0].replace("torsions_", "").replace("torsions", "rotatable")
        ax.set_title(f"t{k}  {'-'.join(atom_names[k])}  ({kind})", fontsize=8.5)
        if k == 0:
            ax.legend(fontsize=7, loc="upper left", frameon=False)
    lengths = ", ".join(f"{label} {s[3]:.0f} ns" for label, s in systems.items())
    fig.suptitle(f"{man['name']} torsions: water vs bound (state 0, all frames; {lengths}; dashed = crystal)",
                 fontsize=12)
    fig.savefig(args.out + ".png", dpi=150)
    plt.close(fig)
    print(f"wrote {args.out}.png and {args.out}.csv")


if __name__ == "__main__":
    main()
