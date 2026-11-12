"""How much of the free (water) ensemble has a structural match in the bound (complex) ensemble, as a
function of the RMSD that counts as a match -- and its area under the curve as one score per ligand.

Clusters: the barrier-cut clusters of the water run and of each complex (report_cut.sbatch, kT_cut 1,
<work>/<c>/<system>-cut/labels.csv), ALL of them (not only those >= 1 %), so the water populations add up to 1.
  d(Y^w, X^p)  = mean symmetry-corrected heavy-atom RMSD between frames of water cluster Y and complex
                 cluster X (as cross_cluster_rmsd.py; at most --per-cluster random frames of each cluster)
  dmin(Y^w)    = min over the complex clusters X of d(Y^w, X^p): how far water cluster Y is from the
                 nearest bound conformer family
  A(r)         = sum of pi^w_Y over the water clusters with dmin(Y^w) < r: the fraction of the free ensemble
                 that has a bound counterpart within r
  score        = (1/r_max) * integral_0^r_max A(r) dr, r_max = --r-max (10 A): 1 = every water conformer
                 matches a bound one at any r (lock and key); lower = more of water is far from anything bound.
Both ensembles in the water ligand's atom order (kl_cross.ensembles). Folate: one curve per complex.

  python -I cluster_match_auc.py --name folic-acid   -> results/<name>/<version>/<name>_cluster_match.png/.csv
  python -I cluster_match_auc.py --plot              -> results/overview/cluster-match/ (all compounds)
"""
import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kl_cross import WORK, REPO, ensembles  # noqa: E402

OUT = os.path.join(REPO, "results", "overview", "cluster-match")
R = np.linspace(0.0, 10.0, 2001)


def labels(cut_dir, n):
    lab = pd.read_csv(f"{cut_dir}/labels.csv", keep_default_na=False, dtype={"cluster": str})["cluster"].to_numpy()
    if len(lab) != n:
        sys.exit(f"{cut_dir}: {len(lab)} labels for {n} frames")
    return lab


def sample(lab, per, rng):
    """{cluster: frame indices}, at most `per` random frames of each cluster."""
    out = {}
    for c in np.unique(lab):
        idx = np.flatnonzero(lab == c)
        out[c] = idx if len(idx) <= per else np.sort(rng.choice(idx, per, replace=False))
    return out


def curve(dmin, mass, r_max):
    acc = np.array([mass[dmin < r].sum() for r in R])
    keep = R <= r_max
    return acc, float(np.trapezoid(acc[keep], R[keep]) / r_max)


def run(args):
    from rdkit import Chem
    from md_tools.analysis import enumerate_symmetry, heavy_atoms, operations_on_atoms, symmetric_rmsd
    man = json.load(open(f"{WORK}/{args.name}/manifest.json"))
    mol = Chem.SDMolSupplier(man["water"]["sdf"], removeHs=False)[0]
    heavy = heavy_atoms(mol)
    kw = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
    ops = operations_on_atoms(enumerate_symmetry(mol, **kw), heavy)
    rng = np.random.default_rng(0)
    rows, curves = [], {}
    for c in man["complexes"]:
        tw, tc, _, _ = ensembles(man, c)
        lw = labels(f"{WORK}/{args.name}/water-cut", tw.n_frames)
        lp = labels(f"{WORK}/{args.name}/{c['label']}-cut", tc.n_frames)
        sw, sp = sample(lw, args.per_cluster, rng), sample(lp, args.per_cluster, rng)
        wn, pn = list(sw), list(sp)
        iw = np.concatenate([sw[k] for k in wn]); ip = np.concatenate([sp[k] for k in pn])
        gw = np.concatenate([[i] * len(sw[k]) for i, k in enumerate(wn)])
        gp = np.concatenate([[j] * len(sp[k]) for j, k in enumerate(pn)])
        D = symmetric_rmsd(tw.xyz[iw][:, heavy] * 10.0, tc.xyz[ip][:, heavy] * 10.0, ops, workers=args.workers)
        Aw = np.zeros((len(wn), len(iw))); Aw[gw, np.arange(len(iw))] = 1
        Ap = np.zeros((len(ip), len(pn))); Ap[np.arange(len(ip)), gp] = 1
        M = (Aw @ D @ Ap) / np.outer(Aw.sum(1), Ap.sum(0))           # cluster x cluster mean RMSD
        mass_w = np.array([(lw == k).mean() for k in wn])
        mass_p = np.array([(lp == k).mean() for k in pn])
        dmin, best = M.min(axis=1), M.argmin(axis=1)
        acc, score = curve(dmin, mass_w, args.r_max)
        curves[c["label"]] = (acc, score)
        for i, k in enumerate(wn):
            rows.append({"complex": c["label"], "water_cluster": k, "water_mass": round(mass_w[i], 6),
                         "dmin_rmsd": round(float(dmin[i]), 4), "best_complex_cluster": pn[best[i]],
                         "best_complex_mass": round(float(mass_p[best[i]]), 6)})
        r50 = R[np.argmax(acc >= 0.5)] if acc.max() >= 0.5 else np.nan
        print(f"{args.name} {c['label']}: {len(wn)} water x {len(pn)} complex clusters, {len(iw)} x {len(ip)} frames; "
              f"score {score:.3f} (r_max {args.r_max:g} A), 50 % of water matched below {r50:.2f} A, "
              f"min dmin {dmin.min():.2f} A", flush=True)
    stem = os.path.join(REPO, "results", args.name, man["version"], f"{args.name}_cluster_match")
    df = pd.DataFrame(rows)
    df.to_csv(stem + ".csv", index=False)
    pd.DataFrame([{"ligand": args.name, "complex": k, "score": round(s, 4), "r_max": args.r_max,
                   "r50": round(float(R[np.argmax(a >= 0.5)]), 3)} for k, (a, s) in curves.items()]).to_csv(
        f"{args.out}/data/cluster_match_{args.name}.csv", index=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5), layout="constrained")
    for (k, (acc, score)), col in zip(curves.items(), ["#b2182b", "#2166ac", "#1b7837", "#762a83"]):
        ax.plot(R, 100 * acc, color=col, lw=2, label=f"{k} complex: score {score:.3f}")
        ax.fill_between(R, 0, 100 * acc, color=col, alpha=0.08)
    ax.set_xlim(0, args.r_max)
    ax.set_ylim(0, 102)
    ax.set_xlabel("RMSD threshold r (Å): a water cluster matches if its nearest complex cluster is closer than r")
    ax.set_ylabel("water population with a bound match (%)")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")
    ax.set_title(f"{args.name}: how much of the free ensemble has a structural counterpart in the bound one\n"
                 f"all barrier-cut clusters (kT_cut 1); mean symmetry-corrected heavy-atom RMSD\n"
                 f"between clusters; score = area / {args.r_max:g} Å", fontsize=9.5)
    fig.savefig(stem + ".png", dpi=150)
    print(f"wrote {stem}.png / .csv")


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from kl_cross import ligand_info
    df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{args.out}/data/cluster_match_*.csv"))], ignore_index=True)
    df = df[~df.ligand.str.endswith("-reduce")]
    info = ligand_info()
    df["rot_bonds"] = [info[l][0] for l in df.ligand]
    df["label"] = [f"{r.ligand} {r.complex}" if (df.ligand == r.ligand).sum() > 1 else r.ligand for r in df.itertuples()]
    df.to_csv(f"{args.out}/data/cluster_match_all.csv", index=False)
    d = df.sort_values("score").reset_index(drop=True)
    lo, hi = int(d.rot_bonds.min()), int(d.rot_bonds.max())
    cmap = plt.get_cmap("viridis", hi - lo + 1)
    norm = BoundaryNorm(np.arange(lo - 0.5, hi + 1.5), cmap.N)
    fig, ax = plt.subplots(figsize=(9, 1.5 + 0.38 * len(d)), layout="constrained")
    ax.barh(range(len(d)), d.score, color=[cmap(norm(v)) for v in d.rot_bonds], edgecolor="k", lw=0.8)
    for i, r in d.iterrows():
        ax.text(r.score + 0.003, i, f"{r.score:.3f}", va="center", fontsize=8)
    ax.set_yticks(range(len(d)), d.label)
    lo_x = max(0.0, float(d.score.min()) - 0.05)
    ax.set_xlim(lo_x, 1.0)
    ax.axvline(1.0, color="0.3", ls=":", lw=1)
    ax.set_xlabel(f"cluster-match score = area under the matched-water-population curve / {df.r_max.iloc[0]:g} Å\n"
                  "(1 = every water conformer has a bound counterpart at any RMSD)")
    ax.grid(alpha=0.3, axis="x")
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, pad=0.02,
                      ticks=range(lo, hi + 1))
    cb.set_label("rotatable bonds")
    ax.set_title("How much of the free ligand's ensemble has a structural match in the bound ensemble\n"
                 "barrier-cut clusters (kT_cut 1), all of them; note the x axis starts at "
                 f"{lo_x:.2f}", fontsize=10)
    fig.savefig(f"{args.out}/cluster_match_scores.png", dpi=150)
    print(f"wrote {args.out}/cluster_match_scores.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--per-cluster", type=int, default=20)
    ap.add_argument("--r-max", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(f"{args.out}/data", exist_ok=True)
    plot(args) if args.plot else run(args)


if __name__ == "__main__":
    main()
