"""Folate RMSD to the crystal conformation for the complex runs WITH crystal waters, and the density
plot against the earlier runs WITHOUT them.

Same measure as FOL-REST2-v4/scripts/rmsd_vs_crystal.py: every frame's folate heavy atoms are
superposed onto the crystal folate (Kabsch) and the lowest RMSD over the symmetry-equivalent atom
assignments is kept (the two O of each carboxylate, the two sides of the benzoyl ring). So it is the
ligand's CONFORMATION, not its position in the pocket (that is check_bound.py's pocket-frame RMSD).
REST2 state 0, every saved frame, heavy atoms only.

New runs (computed here): REST2-complex-xtalwater/FOL-4LRH and FOL-3QL3, folate cut out of the
state-0 solute trajectory. They use the unchanged v4 folate package, so the atom order is the one of
the earlier v4 runs (checked on the elements against packages/FOL/FOL.sdf).
Earlier runs (read, not recomputed; their trajectories are not on sciCORE): the per-frame values in
results/symfirst-2026-10-15/1_v4-new-rest2/folic-acid-v4_rmsd_vs_crystal.csv (v4 water vs both
crystals, v4 4LRH and 3QL3 complexes without crystal waters).

  python -I rmsd_vs_crystal.py [--thesis DIR] [--repo DIR] --out <prefix>
writes <prefix>.csv (new runs, columns run,crystal,time_ns,rmsd_A as the earlier file),
<prefix>.png (new runs over time) and <prefix>_hist.png (density, all six curves).
"""
import argparse
import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from rdkit import Chem

SCICORE = "/scicore/home/lill0000/popa0002/thesis"
OLD_CSV = "results/symfirst-2026-10-15/1_v4-new-rest2/folic-acid-v4_rmsd_vs_crystal.csv"
NEW = {"4LRH": "4LRH complex, crystal waters", "3QL3": "3QL3 complex, crystal waters"}
# (run label, crystal, colour, line style, legend); 4LRH red family, 3QL3 blue family
CURVES = [
    ("water", "4LRH", "#f4a3a3", "-", "folic acid in water vs 4LRH crystal"),
    ("water", "3QL3", "#9ecae1", "-", "folic acid in water vs 3QL3 crystal"),
    ("4LRH complex", "4LRH", "#d62728", "--", "4LRH complex, no crystal waters"),
    ("3QL3 complex", "3QL3", "#1f77b4", "--", "3QL3 complex, no crystal waters"),
    (NEW["4LRH"], "4LRH", "#7f0000", "-", "4LRH complex, crystal waters"),
    (NEW["3QL3"], "3QL3", "#08306b", "-", "3QL3 complex, crystal waters"),
]


def skeleton(mol):
    """Heavy-atom graph, every bond single, no aromaticity or charge: element + connectivity."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rw, new = Chem.RWMol(), {}
    for i in heavy:
        atom = Chem.Atom(mol.GetAtomWithIdx(i).GetAtomicNum())
        atom.SetNoImplicit(True)
        new[i] = rw.AddAtom(atom)
    for b in mol.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in new and c in new:
            rw.AddBond(new[a], new[c], Chem.BondType.SINGLE)
    g = rw.GetMol()
    g.UpdatePropertyCache(strict=False)
    return g, heavy


def kabsch_rmsd(X, Y):
    """RMSD after optimal superposition of every frame of X (F, n, 3) onto Y (n, 3)."""
    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=0)
    H = np.einsum("fni,nj->fij", Xc, Yc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    S[:, -1] *= d                                     # no reflections
    e0 = (Xc ** 2).sum(axis=(1, 2)) + (Yc ** 2).sum()
    return np.sqrt(np.maximum(e0 - 2 * S.sum(axis=1), 0.0) / X.shape[1])


def crystal_orders(sim_g, sim_heavy, path):
    """Crystal heavy-atom coordinates (nm) in simulation order, one array per graph match."""
    ref = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
    ref_g, ref_heavy = skeleton(ref)
    xyz = ref.GetConformer().GetPositions() / 10.0
    matches = sim_g.GetSubstructMatches(ref_g, uniquify=False, useChirality=False, maxMatches=1000)
    assert matches, path
    out = []
    for m in matches:                                  # m[ref_graph_idx] = sim_graph_idx
        inv = {s: r for r, s in enumerate(m)}
        out.append(np.array([xyz[ref_heavy[inv[k]]] for k in range(len(sim_heavy))]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--thesis", default=SCICORE, help="folder holding REST2-complex-xtalwater/")
    ap.add_argument("--repo", default=None, help="bsc-thesis checkout (default: <thesis>/bsc-thesis)")
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()
    repo = args.repo or os.path.join(args.thesis, "bsc-thesis")
    xw = os.path.join(args.thesis, "REST2-complex-xtalwater")

    sim = Chem.MolFromMolFile(os.path.join(xw, "packages/FOL/FOL.sdf"), removeHs=False)
    sim_g, sim_heavy = skeleton(sim)
    elements = [a.GetSymbol().upper() for a in sim.GetAtoms()]

    data = defaultdict(list)                           # (run, crystal) -> (ns, rmsd)
    with open(os.path.join(repo, OLD_CSV)) as f:
        for row in csv.DictReader(f):
            data[(row["run"], row["crystal"])].append((float(row["time_ns"]), float(row["rmsd_A"])))
    data = {k: (np.array([t for t, _ in v]), np.array([r for _, r in v])) for k, v in data.items()}

    rows = []
    for pdb, label in NEW.items():
        run = os.path.join(xw, f"FOL-{pdb}")
        top = md.load(os.path.join(run, "build/built.solute.pdb")).topology
        fol = top.select("resname FOL")
        got = [top.atom(i).element.symbol.upper() for i in fol]
        assert got == elements, f"{pdb}: folate atoms in the topology do not match packages/FOL/FOL.sdf"
        t = md.load(os.path.join(run, "REST2-run1/solute_state0_prod1.nc"),
                    top=os.path.join(run, "build/built.solute.pdb"), atom_indices=fol)
        X = t.xyz[:, sim_heavy, :]
        orders = crystal_orders(sim_g, sim_heavy, os.path.join(repo, f"crystal/fol_{pdb}_A.sdf"))
        r = np.min([kabsch_rmsd(X, Y) for Y in orders], axis=0) * 10.0      # nm -> A
        data[(label, pdb)] = (t.time / 1000.0, r)
        print(f"{label}: {t.n_frames} frames, {len(orders)} symmetry-equivalent assignments tried per frame")

    with open(f"{args.out}.csv", "w") as f:
        f.write("run,crystal,time_ns,rmsd_A\n")
        for pdb, label in NEW.items():
            for a, b in zip(*data[(label, pdb)]):
                f.write(f"{label},{pdb},{a:.3f},{b:.4f}\n")

    # over time: the new runs, with the earlier complexes for reference
    fig, ax = plt.subplots(figsize=(11, 4.8), layout="constrained")
    for run, pdb, colour, style, legend in CURVES[2:]:
        ns, r = data[(run, pdb)]
        ax.scatter(ns, r, s=1.5, color=colour, alpha=0.35, linewidths=0, rasterized=True,
                   marker="s" if style == "--" else "o", label=legend)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("heavy-atom RMSD to crystal conformation (Å)")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    leg = ax.legend(loc="upper right", fontsize=9, markerscale=6)
    for h in leg.legend_handles:
        h.set_alpha(1)
    ax.set_title("Folic acid (v4) in its complexes, with and without crystal waters: RMSD of the ligand "
                 "conformation to the crystal over time\nREST2 state 0; folate superposed on folate, best over "
                 "symmetry-equivalent atoms; one point per saved frame", fontsize=10)
    fig.savefig(f"{args.out}.png", dpi=150)
    plt.close(fig)

    # density: all six curves, same bins, each normalised to area 1
    top_val = max(r.max() for _, r in data.values())
    edges = np.arange(0.0, np.ceil(top_val / 0.05) * 0.05 + 0.05, 0.05)
    fig, ax = plt.subplots(figsize=(9.5, 5.2), layout="constrained")
    for run, pdb, colour, style, legend in CURVES:
        _, r = data[(run, pdb)]
        ax.hist(r, bins=edges, density=True, histtype="step", lw=2.0, color=colour, ls=style,
                label=f"{legend} (n = {len(r)}, median {np.median(r):.2f} Å)")
        rows.append((legend, pdb, len(r), np.median(r), np.mean(r < 1.0), np.mean(r < 1.5)))
    ax.set_xlabel("heavy-atom RMSD to crystal conformation (Å)")
    ax.set_ylabel("probability density (1/Å)")
    ax.set_xlim(0, edges[-1])
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.set_title("Folic acid (v4): distribution of the ligand-conformation RMSD to the crystal\n"
                 "dashed: complexes without crystal waters; solid dark: with crystal waters; light: in water. "
                 "State 0, bins 0.05 Å, area 1", fontsize=10)
    fig.savefig(f"{args.out}_hist.png", dpi=150)
    plt.close(fig)

    print("\n| run | vs crystal | frames | median RMSD (Å) | < 1 Å | < 1.5 Å |")
    print("|---|---|---|---|---|---|")
    for legend, pdb, n, med, f1, f15 in rows:
        print(f"| {legend} | {pdb} | {n} | {med:.2f} | {100 * f1:.0f} % | {100 * f15:.0f} % |")
    print(f"\nwritten {args.out}.csv, {args.out}.png, {args.out}_hist.png")


if __name__ == "__main__":
    main()
