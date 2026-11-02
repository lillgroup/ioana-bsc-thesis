"""Ligand conformation RMSD to the crystal, water vs complex: one density plot per compound.

For every complex in the manifest (prepare.py), two curves against that complex's crystal ligand:
the ligand in water (light) and in the complex (dark). Every frame's heavy atoms are superposed on the
crystal heavy atoms (Kabsch) before the RMSD -- the ligand's SHAPE, not its place in the pocket (that
is binding-check/check_bound.py) -- and the lowest RMSD over the symmetry-equivalent assignments of
the heavy-atom graph is kept (carboxylate O, ring flips). State 0, every saved frame.

  python -I rmsd_density.py --manifest <work>/<compound>/manifest.json --out <prefix>
writes <prefix>.png and <prefix>.csv (median, fraction of frames < 1 and < 1.5 A per curve).
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
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare import skeleton

PALETTE = [("#f4a3a3", "#b2182b"), ("#9ecae1", "#08306b"), ("#a1d99b", "#00441b"), ("#fdd0a2", "#7f2704")]


def kabsch_rmsd(X, Y):
    """RMSD after optimal superposition of every frame of X (F, n, 3) onto Y (n, 3)."""
    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=0)
    H = np.einsum("fni,nj->fij", Xc, Yc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    S[:, -1] *= d
    e0 = (Xc ** 2).sum(axis=(1, 2)) + (Yc ** 2).sum()
    return np.sqrt(np.maximum(e0 - 2 * S.sum(axis=1), 0.0) / X.shape[1])


def rmsd_to_crystal(system, ref):
    """Per-frame heavy-atom RMSD (A) of a system's ligand to one crystal SDF, best over graph matches."""
    mol = Chem.MolFromMolFile(system["sdf"], removeHs=False, sanitize=False)
    g, heavy = skeleton(mol)
    crystal = Chem.MolFromMolFile(ref, removeHs=False, sanitize=False)
    rg, rheavy = skeleton(crystal)
    xyz = crystal.GetConformer().GetPositions() / 10.0
    matches = g.GetSubstructMatches(rg, uniquify=False, useChirality=False, maxMatches=2000)
    if not matches:
        sys.exit(f"{ref}: crystal ligand does not match {system['sdf']}")
    orders = []
    for m in matches:                       # m[crystal graph idx] = system graph idx
        inv = {s: r for r, s in enumerate(m)}
        orders.append(np.array([xyz[rheavy[inv[k]]] for k in range(len(heavy))]))
    t = md.load(system["traj"], top=system["top"])
    X = t.xyz[:, heavy, :]
    return np.min([kabsch_rmsd(X, Y) for Y in orders], axis=0) * 10.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    man = json.load(open(args.manifest))
    if not man["complexes"]:
        sys.exit("no complex: nothing to compare")

    curves = []
    for k, c in enumerate(man["complexes"]):
        light, dark = PALETTE[k % len(PALETTE)]
        ref = c["refs"][0]
        curves.append((f"in water vs {c['label']} crystal", rmsd_to_crystal(man["water"], ref), light, "-"))
        curves.append((f"{c['label']} complex vs {c['label']} crystal", rmsd_to_crystal(c, ref), dark, "-"))

    top = max(r.max() for _, r, _, _ in curves)
    edges = np.arange(0.0, np.ceil(top / 0.05) * 0.05 + 0.05, 0.05)
    fig, ax = plt.subplots(figsize=(9.5, 5.2), layout="constrained")
    rows = []
    for label, r, colour, style in curves:
        ax.hist(r, bins=edges, density=True, histtype="step", lw=2.0, color=colour, ls=style,
                label=f"{label} (n = {len(r)}, median {np.median(r):.2f} Å)")
        rows.append((label, len(r), np.median(r), np.mean(r < 1.0), np.mean(r < 1.5)))
    ax.set_xlabel("heavy-atom RMSD to the crystal conformation (Å)")
    ax.set_ylabel("probability density (1/Å)")
    ax.set_xlim(0, edges[-1])
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_title(f"{man['name']}: ligand conformation RMSD to the crystal, water vs complex\n"
                 "REST2 state 0, ligand superposed on the crystal ligand, best over symmetry-equivalent "
                 "atoms; bins 0.05 Å, area 1", fontsize=10)
    fig.savefig(f"{args.out}.png", dpi=150)
    plt.close(fig)
    with open(f"{args.out}.csv", "w") as f:
        f.write("curve,frames,median_A,frac_below_1A,frac_below_1.5A\n")
        for label, n, med, f1, f15 in rows:
            f.write(f"{label},{n},{med:.3f},{f1:.3f},{f15:.3f}\n")
            print(f"{label:40s} {n:6d} frames  median {med:.2f} A  <1 A {100 * f1:4.0f} %  <1.5 A {100 * f15:4.0f} %")
    print(f"wrote {args.out}.png and {args.out}.csv")


if __name__ == "__main__":
    main()
