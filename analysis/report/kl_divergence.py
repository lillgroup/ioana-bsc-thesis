"""KL divergence between the bound and the free ligand ensemble, on ONE shared clustering:

  KL(pi^p || pi^w) = sum_c pi^p_c ln(pi^p_c / pi^w_c)   (k_B)

pi^p: cluster populations of the ligand in its protein (complex), pi^w: in water. Separate clusterings
of water and complex have unrelated cluster labels, so for every (ligand, complex) the two state-0
ensembles are clustered TOGETHER -- the same number of frames from each, evenly spaced, so neither
decides the cluster boundaries -- with one method, and each environment's populations are counted
on those shared clusters.

  --method cut      Shu-Yu's torsion partitioning (cut_clustering, MD-tools 306b109, env analysis-env-cut),
                    defaults, all frames of the shorter run (and as many, evenly spaced, of the longer)
  --method hdbscan  symmetry-first t_hdbscan, allow_single_cluster=False, MD-tools c30c1c0 (env
                    analysis-env-main), 5000 + 5000 frames (the exact distance matrix is capped at
                    10 000); noise / ambiguous frames are left out and populations renormalised

Reading it. If the complex simply SELECTS one conformer that has population pi in water,
KL = ln(1/pi): selecting a typical water conformer gives KL ~ S_w (the water entropy over the same
clusters), selecting a rare one KL >> S_w, and a conformer water never visits makes KL infinite
(induced fit). So every row also has S_w and the fraction of the complex in clusters with NO water
frame ("unseen in water"). KL uses a pseudocount (--pseudocount frames per cluster, default 0.5) so it
stays finite; the unseen fraction is the direct induced-fit measure. Error bars: block SE over
--blocks consecutive time blocks of each environment (the clustering is not refitted).

The torsions are each system's files from prepare.py (one numbering). The complex trajectory is put
into the water ligand's atom order by the same graph mapping; a tautomeric H with no counterpart (3QL3,
folate 3H) takes the slot of the water's, which no torsion uses.

  python -I kl_divergence.py --method cut [--compounds ...]          -> <out>/kl_<method>.csv
  python -I kl_divergence.py --plot                                   -> <out>/kl_divergence_{cut,hdbscan}.png
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "shuyu_clustering"))
WORK = "/scicore/home/lill0000/popa0002/thesis/analysis-report"
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(REPO, "results", "overview", "kl-divergence")
DEFAULT = ["4j7d", "4ibi", "4ibk", "5mwh", "5ov8", "folic-acid"]       # paracetamol has no complex
KINDS = ("torsions", "torsions_rigid", "torsions_polarH")


def populations(labels, clusters, pseudo=0.0):
    lab = labels[np.isin(labels, clusters)]
    counts = np.array([(lab == c).sum() for c in clusters], float)
    return counts, (counts + pseudo) / (counts + pseudo).sum()


def kl(p, q):
    m = p > 0
    return float((p[m] * np.log(p[m] / q[m])).sum())


def entropy(p):
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def pooled(man, c, method):
    """Torsions (frames x K), joined trajectory (water atom order), environment flag, time order."""
    import mdtraj as md
    from rdkit import Chem
    from cluster_prune import load_torsions
    from prepare import atom_map
    w, cx = man["water"], c
    tw = md.load(w["traj"], top=w["top"])
    tc = md.load(cx["traj"], top=cx["top"])
    n = min(tw.n_frames, tc.n_frames, 5000 if method == "hdbscan" else 10 ** 9)
    iw = np.linspace(0, tw.n_frames - 1, n).round().astype(int)
    ic = np.linspace(0, tc.n_frames - 1, n).round().astype(int)
    qw = list(load_torsions([f"{w['torsions']}/{k}.txt" for k in KINDS if os.path.getsize(f"{w['torsions']}/{k}.txt")]).items())
    qc = list(load_torsions([f"{cx['torsions']}/{k}.txt" for k in KINDS if os.path.getsize(f"{cx['torsions']}/{k}.txt")]).values())
    names, quads = [q[0] for q in qw], [q[1] for q in qw]
    tors = np.vstack([md.compute_dihedrals(tw[iw], quads), md.compute_dihedrals(tc[ic], qc)])
    # complex coordinates in the water ligand's atom order
    amap, _ = atom_map(Chem.MolFromMolFile(w["sdf"], removeHs=False), Chem.MolFromMolFile(cx["sdf"], removeHs=False))
    spare = iter(sorted(set(range(tc.n_atoms)) - {v for v in amap.values() if v is not None}))
    order = [amap[i] if amap.get(i) is not None else next(spare) for i in range(tw.n_atoms)]
    xyz = np.concatenate([tw.xyz[iw], tc.xyz[ic][:, order]])
    joined = md.Trajectory(xyz, tw.topology, time=np.arange(2 * n, dtype=float),
                           unitcell_lengths=None, unitcell_angles=None)
    env = np.array([0] * n + [1] * n)
    return tors, joined, env, names, quads, n


def run_method(args):
    from rdkit import Chem
    import md_tools
    from md_tools.analysis import t_hdbscan, TorsionDefinitions, verify_atom_mapping
    rows = []
    for name in args.compounds:
        man = json.load(open(f"{WORK}/{name}/manifest.json"))
        mol = Chem.SDMolSupplier(man["water"]["sdf"], removeHs=False)[0]
        extra = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
        for c in man["complexes"]:
            tors, joined, env, names, quads, n = pooled(man, c, args.method)
            mapping, why = verify_atom_mapping(mol, [(a.index, b.index) for a, b in joined.topology.bonds],
                                               [a.element.symbol for a in joined.topology.atoms])
            defs = TorsionDefinitions(quads=tuple(quads), names=tuple(names), indexing="zero-based", units="radian")
            fit = t_hdbscan(tors, molecule=mol, torsion_definitions=defs, coordinates=joined, atom_mapping=mapping,
                            **({"allow_single_cluster": False} if args.method == "hdbscan" else {}), **extra)
            if args.method == "cut":
                from md_tools.analysis import cut_clustering
                r = cut_clustering(fit)
                labels, clusters = r["labels"], list(range(r["n_clusters"]))
                settings = {k: v for k, v in r["settings"].items() if k not in ("floor_frames", "reference")}
            else:
                fit = fit.fit()
                labels, clusters = fit.labels_, list(range(fit.n_clusters_))
                s = fit.summary()
                settings = {"route": s["route"], "allow_single_cluster": s["allow_single_cluster"],
                            "mass_floor": s["mass_floor"]}
            lw, lc = labels[env == 0], labels[env == 1]
            cw, pw = populations(lw, clusters, args.pseudocount)
            cc, pc = populations(lc, clusters, args.pseudocount)
            _, pw0 = populations(lw, clusters)
            _, pc0 = populations(lc, clusters)
            unseen = float(pc0[cw == 0].sum())
            # block SE: populations per time block of each environment, same clusters
            vals = []
            for bw, bc in zip(np.array_split(lw, args.blocks), np.array_split(lc, args.blocks)):
                vals.append(kl(populations(bc, clusters, args.pseudocount)[1], populations(bw, clusters, args.pseudocount)[1]))
            row = {"method": args.method, "ligand": name, "complex": c["label"], "frames_each": n,
                   "n_clusters": len(clusters), "KL_kB": round(kl(pc, pw), 4),
                   "KL_block_se": round(float(np.std(vals, ddof=1) / np.sqrt(len(vals))), 4),
                   "S_water_kB": round(entropy(pw0), 4), "S_complex_kB": round(entropy(pc0), 4),
                   "complex_in_unseen_clusters": round(unseen, 4),
                   "water_unassigned": round(float(np.mean(~np.isin(lw, clusters))), 4),
                   "complex_unassigned": round(float(np.mean(~np.isin(lc, clusters))), 4),
                   "top_complex_cluster": int(np.argmax(pc0)),
                   "top_complex_mass": round(float(pc0.max()), 4),
                   "water_mass_of_that_cluster": round(float(pw0[np.argmax(pc0)]), 4),
                   "md_tools": md_tools.__version__, "settings": json.dumps(settings, sort_keys=True, default=str),
                   "pseudocount": args.pseudocount}
            rows.append(row)
            print(f"{name:11s} {c['label']}: {len(clusters):3d} shared clusters on {n}+{n} frames | KL {row['KL_kB']:.2f} "
                  f"+- {row['KL_block_se']:.2f} kB | S_w {row['S_water_kB']:.2f} | complex in clusters unseen in water "
                  f"{100 * unseen:.1f} % | top complex cluster {100 * row['top_complex_mass']:.0f} %, "
                  f"in water {100 * row['water_mass_of_that_cluster']:.1f} %")
    df = pd.DataFrame(rows)
    if len(set(df.settings)) != 1:
        sys.exit(f"settings differ between compounds: {set(df.settings)}")
    os.makedirs(args.out, exist_ok=True)
    df.to_csv(f"{args.out}/kl_{args.method}.csv", index=False)
    print(f"wrote {args.out}/kl_{args.method}.csv")


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    titles = {"cut": "torsion partitioning (barrier cuts)", "hdbscan": "symmetry-first HDBSCAN (no single cluster)"}
    for method, title in titles.items():
        f = f"{args.out}/kl_{method}.csv"
        if not os.path.exists(f):
            print(f"skip {method}: no {f}")
            continue
        df = pd.read_csv(f)
        df["label"] = df.ligand + "\n" + df.complex
        x = np.arange(len(df))
        fig, ax = plt.subplots(figsize=(11.5, 5.8), layout="constrained")
        ax.bar(x - 0.18, df.KL_kB, 0.36, yerr=df.KL_block_se, capsize=3, color="#b2182b",
               label=r"KL($\pi^{complex}\,\|\,\pi^{water}$)")
        ax.bar(x + 0.18, df.S_water_kB, 0.36, color="#9a9a9a",
               label=r"$S_w$ = KL of selecting a typical water conformer")
        for xi, (_, r) in zip(x, df.iterrows()):
            ax.annotate(f"{100 * r.complex_in_unseen_clusters:.0f} %\nunseen", (xi - 0.18, r.KL_kB + r.KL_block_se),
                        textcoords="offset points", xytext=(0, 4), ha="center", fontsize=7.5, color="#b2182b")
            ax.annotate(f"top: {100 * r.top_complex_mass:.0f} % bound,\n{100 * r.water_mass_of_that_cluster:.1f} % in water",
                        (xi, 0), textcoords="offset points", xytext=(0, -38), ha="center", fontsize=7, color="0.3",
                        annotation_clip=False)
        ax.set_xticks(x, df.label, fontsize=9)
        ax.tick_params(axis="x", pad=40)
        ax.set_ylabel("k$_B$")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(loc="upper left", fontsize=9)
        ax.set_title(f"Bound vs free ligand ensemble: KL divergence on one shared clustering of water + complex frames\n"
                     f"{title}; KL ~ S_w: the protein selects a conformer water samples (conformational selection); "
                     f"KL >> S_w or\nmass in clusters unseen in water: it binds conformers water rarely or never visits "
                     f"(induced fit). Pseudocount {df.pseudocount.iloc[0]:g} frames/cluster; bars: block SE",
                     fontsize=9.5)
        fig.savefig(f"{args.out}/kl_divergence_{method}.png", dpi=150)
        plt.close(fig)
        print(f"wrote {args.out}/kl_divergence_{method}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--method", choices=("cut", "hdbscan"))
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--compounds", nargs="+", default=DEFAULT)
    ap.add_argument("--pseudocount", type=float, default=0.5)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    if args.plot:
        plot(args)
    elif args.method:
        run_method(args)
    else:
        ap.error("--method or --plot")


if __name__ == "__main__":
    main()
