"""Bioactive-conformer analysis: how close does each ensemble, and each of its clusters, come to the
crystal conformation of the ligand?

For ONE crystal ligand and one or more ensembles (e.g. water and protein), every frame's ligand heavy
atoms are superposed onto the crystal ligand (Kabsch) and the lowest RMSD over the
symmetry-equivalent atom assignments is kept (carboxylate O, ring sides, CF3 F, ...). So it measures
the ligand's CONFORMATION, not its position in the pocket. Per ensemble:
  RMSD-<ens>       the whole ensemble (every frame of the trajectory)
  RMSD-<ens>-<c>   the frames cluster_prune.py put in cluster <c> (its labels.csv, matched on time;
                   noise / ambiguous / unassigned frames form their own "unassigned" curve)

  python -I bioactive_rmsd.py --sdf FOL.sdf --top built.solute.pdb --crystal fol_4LRH_A.sdf \
      --ensemble water   water/solute_state0_prod1.nc  <water clustering>/labels.csv \
      --ensemble protein 4LRH/FOL_state0_prod1.nc      <complex clustering>/labels.csv \
      --out results/bioactive/folic-acid-v4_4LRH
(LABELS "-" = no clustering, whole ensemble only.) The trajectory must contain the --sdf ligand's
atoms in --sdf order (state-0 solute or ligand-only trajectory; --select picks the ligand otherwise).
Writes <out>.csv (per frame: ensemble, time_ps, cluster, rmsd_A), <out>_summary.csv (per ensemble
and cluster: frames, population, min / median RMSD, fraction below 1 / 1.5 / 2 A) and <out>.png
(top: the ensemble density with each cluster's share, area = population; bottom: each cluster's
own density, area 1).
"""
import argparse
import csv
import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from rdkit import Chem

# categorical slots in fixed order (validated default palette); the 8th is "other clusters"
CLUSTER_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"]
OTHER_COLOUR, UNASSIGNED_COLOUR, ENSEMBLE_COLOUR = "#e34948", "#8c8c8c", "#3a3a3a"
UNASSIGNED = {"noise", "ambiguous", "unassigned"}
THRESHOLDS = (1.0, 1.5, 2.0)


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


def crystal_orders(sim_g, sim_heavy, path, max_matches):
    """Crystal heavy-atom coordinates (nm) in simulation order, one array per distinct graph match."""
    ref = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
    ref_g, ref_heavy = skeleton(ref)
    xyz = ref.GetConformer().GetPositions() / 10.0
    matches = sim_g.GetSubstructMatches(ref_g, uniquify=False, useChirality=False, maxMatches=max_matches)
    assert matches, f"{path}: crystal ligand does not match the --sdf ligand (heavy-atom graph)"
    if len(matches) == max_matches:
        raise SystemExit(f"{path}: {max_matches} graph matches, search may be incomplete; raise --max-matches")
    out = []
    for m in matches:                                  # m[ref_graph_idx] = sim_graph_idx
        inv = {s: r for r, s in enumerate(m)}
        out.append(np.array([xyz[ref_heavy[inv[k]]] for k in range(len(sim_heavy))]))
    return out


def best_rmsd(X, orders, chunk=2000):
    """Lowest RMSD (A) over the crystal atom assignments, per frame, in frame chunks."""
    out = np.empty(len(X))
    for s in range(0, len(X), chunk):
        out[s:s + chunk] = np.min([kabsch_rmsd(X[s:s + chunk], Y) for Y in orders], axis=0) * 10.0
    return out


def read_labels(path, time_ps):
    """Cluster name per trajectory frame from a labels.csv (matched on time; '' = frame not clustered)."""
    lab = np.full(len(time_ps), "", dtype=object)
    if path == "-":
        return lab
    index = {round(float(t), 1): i for i, t in enumerate(time_ps)}
    missing = 0
    with open(path) as f:
        for row in csv.DictReader(f):
            i = index.get(round(float(row["time_ps"]), 1))
            if i is None:
                missing += 1
                continue
            lab[i] = "unassigned" if row["label"] in UNASSIGNED else row["label"]
    if missing:
        raise SystemExit(f"{path}: {missing} labelled times not in the trajectory (wrong trajectory?)")
    return lab


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sdf", required=True, help="ligand as simulated (atom order of the trajectory)")
    ap.add_argument("--top", required=True, help="topology of the trajectories (e.g. built.solute.pdb)")
    ap.add_argument("--select", default=None,
                    help="mdtraj selection of the ligand in --top (default: all atoms)")
    ap.add_argument("--crystal", required=True, help="crystal ligand (sdf, with coordinates)")
    ap.add_argument("--crystal-name", default=None, help="label in titles (default: file name)")
    ap.add_argument("--ensemble", nargs=3, action="append", required=True, metavar=("NAME", "TRAJ", "LABELS"),
                    help="ensemble name, its trajectory, its cluster_prune.py labels.csv (or '-'); repeat")
    ap.add_argument("--name", default="", help="system name for the titles")
    ap.add_argument("--max-matches", type=int, default=100000, help="cap on the symmetry assignments")
    ap.add_argument("--bin", type=float, default=0.05, help="histogram bin width (A)")
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()
    cname = args.crystal_name or os.path.splitext(os.path.basename(args.crystal))[0]

    sim = Chem.MolFromMolFile(args.sdf, removeHs=False)
    sim_g, sim_heavy = skeleton(sim)
    orders = crystal_orders(sim_g, sim_heavy, args.crystal, args.max_matches)
    top = md.load_topology(args.top)
    sel = top.select(args.select) if args.select else np.arange(top.n_atoms)
    got = [top.atom(i).element.symbol.upper() for i in sel]
    want = [a.GetSymbol().upper() for a in sim.GetAtoms()]
    assert got == want, "ligand atoms in --top (after --select) do not match --sdf in order and element"
    print(f"{cname}: {len(orders)} symmetry-equivalent atom assignments tried per frame")

    ens = []                                           # (name, time_ps, rmsd, labels)
    for name, traj, labels in args.ensemble:
        t = md.load(traj, top=args.top, atom_indices=sel)
        r = best_rmsd(t.xyz[:, sim_heavy, :], orders)
        lab = read_labels(labels, t.time)
        ens.append((name, t.time, r, lab))
        print(f"{name}: {t.n_frames} frames, {np.sum(lab != '')} labelled, median RMSD {np.median(r):.2f} A")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(f"{args.out}.csv", "w") as f:
        f.write("ensemble,time_ps,cluster,rmsd_A\n")
        for name, time, r, lab in ens:
            for a, c, b in zip(time, lab, r):
                f.write(f"{name},{a:.1f},{c},{b:.4f}\n")

    def stats(name, key, r, pop):
        return dict(ensemble=name, cluster=key, frames=len(r), population=pop, min_A=r.min(),
                    median_A=np.median(r), **{f"below_{x:g}A": np.mean(r < x) for x in THRESHOLDS})

    rows, panels = [], []
    for name, time, r, lab in ens:
        rows.append(stats(name, "all", r, 1.0))
        done = lab != ""
        counts = Counter(lab[done])
        clusters = sorted(c for c in counts if c != "unassigned")
        for c in clusters:
            m = done & (lab == c)
            rows.append(stats(name, c, r[m], m.sum() / done.sum()))
        # at most 7 drawn on their own: the crystal-like ones (median < 1 A) first, then by population;
        # the rest fold into "other". Colours by name among the drawn ones (A is slot 1).
        crystal_like = lambda c: np.median(r[done & (lab == c)]) < 1.0
        shown = sorted(clusters, key=lambda c: (not crystal_like(c), -counts[c]))[:len(CLUSTER_COLOURS)]
        rest = [c for c in clusters if c not in shown]
        curves = []                                    # (label, colour, rmsd values, population)
        for c, colour in zip(sorted(shown), CLUSTER_COLOURS):
            m = done & (lab == c)
            curves.append((c, colour, r[m], m.sum() / done.sum()))
        if rest:
            m = done & np.isin(lab, rest)
            curves.append((f"other ({', '.join(rest)})", OTHER_COLOUR, r[m], m.sum() / done.sum()))
        if counts.get("unassigned"):
            m = done & (lab == "unassigned")
            rows.append(stats(name, "unassigned", r[m], m.sum() / done.sum()))
            curves.append(("unassigned", UNASSIGNED_COLOUR, r[m], m.sum() / done.sum()))
        panels.append((name, r, curves))

    with open(f"{args.out}_summary.csv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for row in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in row.items()})

    top_val = max(r.max() for _, r, _ in panels)
    edges = np.arange(0.0, np.ceil(top_val / args.bin) * args.bin + args.bin, args.bin)
    centres = 0.5 * (edges[1:] + edges[:-1])
    fig, axes = plt.subplots(2, len(panels), figsize=(6.2 * len(panels), 8.4), sharex=True,
                             squeeze=False, layout="constrained")
    for j, (name, r, curves) in enumerate(panels):
        ax = axes[0, j]
        dens, _ = np.histogram(r, bins=edges, density=True)
        ax.fill_between(centres, dens, step="mid", color="#d9d9d9", lw=0)
        ax.step(centres, dens, where="mid", color=ENSEMBLE_COLOUR, lw=2,
                label=f"RMSD-{name} (all {len(r)} frames, median {np.median(r):.2f} Å)")
        for c, colour, rc, pop in curves:              # share of the ensemble: area = population
            h, _ = np.histogram(rc, bins=edges, density=True)
            ax.step(centres, h * pop, where="mid", color=colour, lw=1.6,
                    ls="--" if c == "unassigned" else "-",
                    label=f"{c}: {100 * pop:.1f} %, median {np.median(rc):.2f} Å")
        ax.set_title(f"RMSD-{name}: whole ensemble and each cluster's share", fontsize=10)
        ax.legend(fontsize=7.5, loc="upper right")
        ax = axes[1, j]
        for c, colour, rc, pop in curves:              # each cluster on its own, area 1
            h, _ = np.histogram(rc, bins=edges, density=True)
            ax.step(centres, h, where="mid", color=colour, lw=1.6, ls="--" if c == "unassigned" else "-",
                    label=f"RMSD-{name}-{c} (n = {len(rc)}, < 1 Å: {100 * np.mean(rc < 1):.0f} %)")
        ax.set_title(f"RMSD-{name}-<cluster>: each cluster normalised to area 1", fontsize=10)
        if curves:
            ax.legend(fontsize=7.5, loc="upper right")
        else:
            ax.text(0.5, 0.5, "no clustering given", transform=ax.transAxes, ha="center", color="#8c8c8c")
        ax.set_xlabel(f"heavy-atom RMSD to the {cname} crystal conformation (Å)")
    for ax in axes.flat:
        ax.set_xlim(0, edges[-1])
        ax.set_ylim(0, ax.get_ylim()[1] * 1.5)          # headroom for the legend
        ax.set_ylabel("probability density (1/Å)")
        ax.grid(alpha=0.3)
    fig.suptitle(f"{args.name} bioactive-conformer analysis vs the {cname} crystal\nligand superposed on "
                 f"the crystal ligand, best over symmetry-equivalent atoms; REST2 state 0, bins {args.bin:g} Å"
                 "\ncluster shares from the clustering's labelled (strided) frames; noise + ambiguous = "
                 "unassigned", fontsize=10)
    fig.savefig(f"{args.out}.png", dpi=150)
    plt.close(fig)

    print(f"\n| ensemble | cluster | frames | population | min (Å) | median (Å) | < 1 Å | < 1.5 Å | < 2 Å |")
    print("|---|---|---|---|---|---|---|---|---|")
    for x in rows:
        print(f"| {x['ensemble']} | {x['cluster']} | {x['frames']} | {100 * x['population']:.1f} % | "
              f"{x['min_A']:.2f} | {x['median_A']:.2f} | " +
              " | ".join(f"{100 * x[f'below_{t:g}A']:.0f} %" for t in THRESHOLDS) + " |")
    print(f"\nwritten {args.out}.csv, {args.out}_summary.csv, {args.out}.png")


if __name__ == "__main__":
    main()
