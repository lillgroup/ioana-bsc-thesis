"""Cost of the bioactive conformation in water (Zivanovic et al., JCTC 16, 6575 (2020), eqs. 2-3):
how rare the crystal (bound) conformation is in the free ligand's REST2 ensemble.

  dG_dist   = -ln P_b         P_b = population of the WATER barrier-cut cluster the crystal falls in
  dG_strain = -ln(P_b / P_0)  P_0 = population of the most populated water cluster (0: the bound
                              conformation is the free ligand's favourite one)
in k_B T. Cluster and crystal placement: report_cut.sbatch's water run (<work>/<c>/water-cut:
summary.json crystal_assignment / cluster_mass, labels.csv), the gold frame of *_water_clusters_cut.png.
Crystals that cluster_cut could not place (a cut through a polar H the crystal has no position for, or
symmetry-equivalent atom assignments in different clusters) are left out of the figures (listed in the
CSV and the figure note). Crystals of one compound in the same water cluster are one row.

Cluster-free cross-check: P_b(r) = fraction of water frames within r = 1.0 / 1.5 / 2.0 A heavy-atom
RMSD of the crystal (superposed, minimised over the molecule's symmetry operations), dG = -ln P_b(r).
Errors: block SE over --blocks contiguous time blocks of the water run, propagated to the log
(sigma_dG = SE(P)/P; strain: both populations, covariance ignored).

Figures (results/overview/bioactive-cost/):
  bioactive_cost_ranked.png      one row per ligand, sorted by dG_dist: dG_dist (dot, block SE),
                                 dG_strain (hollow diamond), RMSD-based dG from 2.0 to 1.0 A (grey bar,
                                 tick at 1.5 A); ligand names coloured by the barrier-cut KL regime
  bioactive_cost_cumulative.png  % of ligands with dG_dist <= x
  bioactive_cost.csv

  python -I bioactive_cost.py [--compounds ...] [--blocks 5]
"""
import argparse
import csv
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
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "shuyu_clustering"))
from cluster_prune import _heavy_graph  # noqa: E402
from md_tools.analysis import (verify_atom_mapping, enumerate_symmetry, heavy_atoms,  # noqa: E402
                               operations_on_atoms, symmetric_rmsd)

WORK = "/scicore/home/lill0000/popa0002/thesis/analysis-report"
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(REPO, "results", "overview", "bioactive-cost")
DEFAULT = [r["name"] for r in csv.DictReader(open(os.path.join(HERE, "compounds.csv")))]
CUTOFFS = (1.0, 1.5, 2.0)
REGIMES = {"lock and key": "#2ca25f", "conformational selection": "#3182bd", "induced fit": "#e6550d"}


def crystal_heavy_xyz(ref_sdf, mol, X0, operations):
    """Crystal heavy-atom coordinates (A) in the simulated ligand's heavy-atom order. The heavy-atom
    graph match is not unique (symmetry, and element-only matching); the match closest to frame 0
    after the symmetry operations is kept -- any symmetry-equivalent one gives the same RMSD."""
    ref = Chem.MolFromMolFile(ref_sdf, removeHs=False, sanitize=False)
    sim_g, sim_heavy = _heavy_graph(mol)
    ref_g, ref_heavy = _heavy_graph(ref)
    matches = sim_g.GetSubstructMatches(ref_g, uniquify=False, useChirality=False, maxMatches=1000)
    if not matches:
        sys.exit(f"{ref_sdf}: heavy-atom graph does not match the simulated ligand")
    xyz = ref.GetConformer().GetPositions()
    cands = []
    for m in matches:                                  # m[ref graph idx] = sim graph idx
        inv = {s: r for r, s in enumerate(m)}
        cands.append(xyz[[ref_heavy[inv[k]] for k in range(len(sim_heavy))]])
    d = symmetric_rmsd(np.asarray(cands), X0[None], operations)[:, 0]
    return cands[int(np.argmin(d))]


def block_pop(mask, blocks):
    """Population of a boolean frame mask and its block SE over contiguous time blocks."""
    parts = np.array_split(np.asarray(mask, float), blocks)
    vals = np.array([p.mean() for p in parts])
    return float(np.mean(mask)), float(vals.std(ddof=1) / np.sqrt(blocks))


def regimes():
    path = os.path.join(REPO, "results", "overview", "kl-divergence", "data", "kl_cross_cut.csv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path)
    return {(r.ligand, r.complex): r.regime for r in df.itertuples()} if "regime" in df else {}


def compute(name, blocks, reg):
    man = json.load(open(f"{WORK}/{name}/manifest.json"))
    wdir = f"{WORK}/{name}/water-cut"
    summary = json.load(open(f"{wdir}/summary.json"))
    labels = pd.read_csv(f"{wdir}/labels.csv", keep_default_na=False, dtype={"cluster": str})["cluster"].to_numpy()
    masses = summary["cluster_mass"]
    top = max(masses, key=masses.get)
    w = man["water"]
    t = md.load(w["traj"], top=w["top"], stride=int(summary["stride"]))
    if t.n_frames != len(labels):
        sys.exit(f"{name}: {t.n_frames} water frames, labels.csv has {len(labels)}")
    mol = Chem.SDMolSupplier(w["sdf"], removeHs=False)[0]
    mapping, _ = verify_atom_mapping(mol, [(a.index, b.index) for a, b in t.topology.bonds],
                                     [a.element.symbol for a in t.topology.atoms])
    heavy = heavy_atoms(mol)
    kw = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
    ops = operations_on_atoms(enumerate_symmetry(mol, **kw), heavy)
    X = t.xyz[:, np.asarray(mapping)[heavy], :] * 10.0
    complexes = [c["label"] for c in man["complexes"]]
    rows, seen = [], set()
    for ref in w["refs"]:
        stem = os.path.splitext(os.path.basename(ref))[0]
        a = summary.get("crystal_assignment", {}).get(stem, {"cluster": None, "note": "no assignment"})
        c = a["cluster"]
        if c is not None and c in seen:
            continue
        seen.add(c)
        cplx = next((x for x in complexes if x.lower() in stem.lower()), complexes[0] if complexes else "")
        rmsd = symmetric_rmsd(X, crystal_heavy_xyz(ref, mol, X[0], ops)[None], ops)[:, 0]
        row = {"ligand": name, "crystal": stem, "complex": cplx, "regime": reg.get((name, cplx), ""),
               "water_cluster": c, "placed": c is not None, "note": a.get("note", ""),
               "n_water_clusters": len(masses), "frames": len(labels), "min_rmsd_A": round(float(rmsd.min()), 3)}
        if c is not None:
            pb, sb = block_pop(labels == c, blocks)
            p0, s0 = block_pop(labels == top, blocks)
            row.update(P_b=pb, P_b_se=sb, top_cluster=top, P_0=p0, dG_dist=-np.log(pb), dG_dist_se=sb / pb,
                       dG_strain=-np.log(pb / p0), dG_strain_se=float(np.hypot(sb / pb, s0 / p0)))
        for r in CUTOFFS:
            p, s = block_pop(rmsd <= r, blocks)
            row[f"P_rmsd_{r:g}A"] = p
            row[f"dG_rmsd_{r:g}A"] = -np.log(p) if p > 0 else np.inf
        rows.append(row)
        print(f"{name:11s} {stem:14s} " + (f"cluster {c:4s} P_b {100 * row['P_b']:6.2f} %  dG_dist {row['dG_dist']:5.2f} "
              f"+- {row['dG_dist_se']:.2f}  dG_strain {row['dG_strain']:5.2f}" if c else f"not placed ({a.get('note', '')})")
              + f" | min RMSD {row['min_rmsd_A']:.2f} A, within 1/1.5/2 A: "
              + "/".join(f"{100 * row[f'P_rmsd_{r:g}A']:.2f}" for r in CUTOFFS) + " %", flush=True)
    return rows


def plot_ranked(df, left_out, path):
    d = df.sort_values("dG_dist").reset_index(drop=True)
    n = len(d)
    xmax = max(8.0, float(np.nanmax(np.r_[d.dG_dist + d.dG_dist_se, d["dG_rmsd_2A"].replace(np.inf, np.nan)])) + 0.8)
    fig, ax = plt.subplots(figsize=(10, 2.2 + 0.42 * n), layout="constrained")
    for x0, x1, c in ((0, 3, "#3182bd"), (3, 5, "#fdae6b"), (5, xmax, "#e6550d")):
        ax.axvspan(x0, x1, color=c, alpha=0.06, lw=0)
    for x in (3, 5):
        ax.axvline(x, color="0.4", ls="--", lw=1)
    for i, r in d.iterrows():
        lo, mid, hi = r["dG_rmsd_2A"], r["dG_rmsd_1.5A"], r["dG_rmsd_1A"]
        if np.isfinite(lo):
            ax.plot([lo, min(hi, xmax)], [i, i], color="0.75", lw=7, solid_capstyle="butt", zorder=1)
            if not np.isfinite(hi):
                ax.annotate("", (xmax, i), (xmax - 0.35, i), arrowprops=dict(arrowstyle="->", color="0.55"))
            if np.isfinite(mid):
                ax.plot([mid], [i], "|", color="0.35", ms=11, mew=1.5, zorder=2)
        ax.plot(r.dG_strain, i, "D", mfc="white", mec="k", ms=7, mew=1.3, zorder=3)
        ax.errorbar(r.dG_dist, i, xerr=r.dG_dist_se, fmt="o", color=REGIMES.get(r.regime, "k"), mec="k",
                    ms=9, capsize=3, ecolor="k", zorder=4)
    ax.set_yticks(range(n), [f"{r.ligand} ({r.complex})" if (d.ligand == r.ligand).sum() > 1 else r.ligand
                             for r in d.itertuples()])
    for tick, r in zip(ax.get_yticklabels(), d.itertuples()):
        tick.set_color(REGIMES.get(r.regime, "k"))
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel(r"$\Delta G$  (k$_B$T)")
    top = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
    pcts = [50, 10, 1, 0.1, 0.01]
    top.set_xticks([-np.log(p / 100) for p in pcts if -np.log(p / 100) <= xmax],
                   [f"{p:g} %" for p in pcts if -np.log(p / 100) <= xmax])
    top.set_xlabel("population in water, e$^{-\\Delta G}$")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([], [], marker="o", ls="", color="0.5", mec="k", ms=8,
                      label=r"$\Delta G_{dist} = -\ln P_b$  (crystal's water cluster; bar: block SE)"),
               Line2D([], [], marker="D", ls="", mfc="white", mec="k", ms=7,
                      label=r"$\Delta G_{strain} = -\ln(P_b/P_0)$  (vs the most populated water cluster)"),
               Line2D([], [], color="0.75", lw=7, label="cluster-free: frames within 2.0 → 1.0 Å RMSD of the crystal (tick: 1.5 Å)")]
    handles += [Patch(color=c, label=f"KL regime: {k}") for k, c in REGIMES.items()]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=8, frameon=False)
    ax.set_title("Cost of the bioactive (crystal) conformation for the free ligand in water\n"
                 "barrier-cut water clusters, REST2 state 0; dashed: 3 and 5 k$_B$T (Zivanovic et al. 2020)"
                 + (f"\nleft out (crystal not placed in a water cluster): {', '.join(left_out)}" if left_out else ""),
                 fontsize=10)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_cumulative(df, left_out, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5), layout="constrained")
    xmax = max(8.0, float(df.dG_dist.max()) + 0.5)
    grid = np.linspace(0, xmax, 801)
    v = df["dG_dist"].to_numpy()
    ax.plot(grid, [100 * np.mean(v <= g) for g in grid], color="k", lw=2,
            label=r"$\Delta G_{dist} = -\ln P_b$  (crystal's water cluster)")
    for x in (3, 5):
        ax.axvline(x, color="0.6", ls="--", lw=1)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 102)
    ax.set_xlabel(r"$\Delta G$  (k$_B$T)")
    ax.set_ylabel(f"% of ligands with ΔG ≤ x  (n = {len(df)})")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.5)
    ax.set_title("Cumulative cost of the bioactive conformation in water\n"
                 "barrier-cut water clusters; dashed: 3 and 5 k$_B$T (Zivanovic et al. 2020)"
                 + (f"\nleft out: {', '.join(left_out)}" if left_out else ""), fontsize=9.5)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--compounds", nargs="+", default=DEFAULT)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    reg = regimes()
    rows = [r for name in args.compounds for r in compute(name, args.blocks, reg)]
    df = pd.DataFrame(rows)
    df.round(5).to_csv(f"{args.out}/bioactive_cost.csv", index=False)
    placed = df[df.placed].copy()
    left_out = sorted(set(df[~df.placed].ligand) - set(placed.ligand))
    plot_ranked(placed, left_out, f"{args.out}/bioactive_cost_ranked.png")
    plot_cumulative(placed, left_out, f"{args.out}/bioactive_cost_cumulative.png")
    print(f"wrote {args.out}/bioactive_cost.csv, bioactive_cost_ranked.png, bioactive_cost_cumulative.png "
          f"({len(placed)} rows; left out: {', '.join(left_out) or 'none'})")


if __name__ == "__main__":
    main()
