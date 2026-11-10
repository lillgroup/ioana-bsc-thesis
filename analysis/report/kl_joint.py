"""KL and Jensen-Shannon divergence between the free (water) and bound (complex) ligand ensembles on ONE
joint clustering: barrier cuts (cut_clustering, MD-tools 306b109, defaults) on the water and complex
frames together, then each ensemble's populations pi^w, pi^p on those clusters.

  equal weight   the longer run is randomly subsampled (seed 0) to the shorter one's frame count for the
                 clustering, so neither ensemble decides the cuts alone; the clustered ensemble is then
                 m = (pi^w + pi^p) / 2. EVERY frame of both runs is labelled afterwards by walking the
                 joint cut tree (kl_cross.cut_assign) on MD-tools' symmetry-invariant coordinates, so the
                 populations use all frames.
  KL             KL(p||w) and KL(w||p) on the joint clusters. A cluster holding frames of only one
                 ensemble makes one direction infinite: populations get a pseudocount of 0.5 frame per
                 cluster ((n_c + 0.5) / (N + 0.5 k)); the CSV gives the number of such clusters and the
                 mass of the other ensemble in them.
  JS             JS = (KL(w||m) + KL(p||m)) / 2 with m = (pi^w + pi^p) / 2 (no pseudocount: m > 0 wherever
                 either has frames); each half lies in [0, ln 2] k_B, ln 2 = no cluster shared.
  errors         block SE over --blocks paired time blocks (block i of water with block i of the complex).
Folate: one joint clustering per complex (water + 4LRH, water + 3QL3). Torsions and atom order: the water
ligand's (kl_cross.ensembles), so t<k> and the symmetry are the same as in every other figure.

  python -I kl_joint.py --compounds 4ibi        -> <out>/kl_joint_cut_4ibi.csv   (env analysis-env-cut)
  python -I kl_joint.py --plot                  -> <out>/kl_joint_cut.png, js_joint_halves_cut.png,
                                                   js_joint_ranked_cut.png, kl_joint_cut.csv (env main)
  bash kl_joint_submit.sh: one 30-min task per compound, then the figures.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from kl_cross import (WORK, OUT, DEFAULT, ensembles, thdbscan, cut_assign,  # noqa: E402
                      ligand_info, point_label)

LN2 = float(np.log(2.0))


def counts(labels, k):
    return np.bincount(labels[(labels >= 0) & (labels < k)], minlength=k).astype(float)


def kl_pc(nq, nr, pc=0.5):
    """KL(q || r) from cluster counts, with a pseudocount of pc frames per cluster on both."""
    q = (nq + pc) / (nq.sum() + pc * len(nq))
    r = (nr + pc) / (nr.sum() + pc * len(nr))
    return float((q * np.log(q / r)).sum())


def js_halves(nw, npx):
    w, p = nw / max(nw.sum(), 1), npx / max(npx.sum(), 1)
    m = 0.5 * (w + p)
    half = lambda a: float((a[a > 0] * np.log(a[a > 0] / m[a > 0])).sum())
    return half(w), half(p)


def measures(nw, npx):
    jw, jp = js_halves(nw, npx)
    return {"KL_p_w": kl_pc(npx, nw), "KL_w_p": kl_pc(nw, npx), "JS_w_m": jw, "JS_p_m": jp, "JS": 0.5 * (jw + jp)}


def run(args):
    import mdtraj as md
    import md_tools
    from rdkit import Chem
    from md_tools.analysis import cut_clustering
    from md_tools.analysis._t_cut import _invariant_coordinates
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(0)
    for name in args.compounds:
        man = json.load(open(f"{WORK}/{name}/manifest.json"))
        mol = Chem.SDMolSupplier(man["water"]["sdf"], removeHs=False)[0]
        extra = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
        rows = []
        for c in man["complexes"]:
            tw, tc, names, quads = ensembles(man, c)
            tw = md.Trajectory(tw.xyz, tw.topology, time=tw.time)    # no box, as the reordered complex (to join)
            n = min(tw.n_frames, tc.n_frames)
            iw = np.sort(rng.choice(tw.n_frames, n, replace=False)) if tw.n_frames > n else np.arange(n)
            ic = np.sort(rng.choice(tc.n_frames, n, replace=False)) if tc.n_frames > n else np.arange(n)
            joint = md.join([tw[iw], tc[ic]])
            r = cut_clustering(thdbscan(joint, mol, names, quads, extra))
            k = r["n_clusters"]
            # every frame of both runs, through the joint tree (one symmetry enumeration for both)
            full = md.join([tw, tc])
            Y, cname = _invariant_coordinates(thdbscan(full, mol, names, quads, extra), r["settings"]["z_min"])[:2]
            if list(cname) != list(r["coordinates"]):
                sys.exit(f"{name} {c['label']}: cut coordinates differ: {cname} vs {r['coordinates']}")
            lab = cut_assign(r, Y)
            lw, lc = lab[:tw.n_frames], lab[tw.n_frames:]
            same = float(np.mean(np.r_[lw[iw], lc[ic]] == r["labels"]))      # tree walk == cut_clustering's labels
            nw, npx = counts(lw, k), counts(lc, k)
            res = measures(nw, npx)
            blocks = [measures(counts(bw, k), counts(bc, k))
                      for bw, bc in zip(np.array_split(lw, args.blocks), np.array_split(lc, args.blocks))]
            se = {key: float(np.std([b[key] for b in blocks], ddof=1) / np.sqrt(args.blocks)) for key in res}
            only_p, only_w = (nw == 0) & (npx > 0), (npx == 0) & (nw > 0)
            row = {"method": "cut", "ligand": name, "complex": c["label"], "n_joint_clusters": int(k),
                   "frames_water": tw.n_frames, "frames_complex": tc.n_frames, "frames_each_in_clustering": int(n),
                   "tree_walk_agreement": round(same, 4),
                   "complex_only_clusters": int(only_p.sum()), "complex_mass_in_them": round(float(npx[only_p].sum() / npx.sum()), 4),
                   "water_only_clusters": int(only_w.sum()), "water_mass_in_them": round(float(nw[only_w].sum() / nw.sum()), 4)}
            for key in res:
                row[key] = round(res[key], 4)
                row[f"{key}_se"] = round(se[key], 4)
            row.update(md_tools=md_tools.__version__,
                       settings=json.dumps({k_: v for k_, v in r["settings"].items() if k_ not in ("floor_frames", "reference")},
                                           sort_keys=True, default=str))
            rows.append(row)
            print(f"{name:11s} {c['label']}: {k:3d} joint clusters (tree walk = labels for {100 * same:.2f} % of the "
                  f"clustered frames) | KL(p||w) {res['KL_p_w']:.2f}  KL(w||p) {res['KL_w_p']:.2f} | JS {res['JS']:.3f} "
                  f"(w {res['JS_w_m']:.3f}, p {res['JS_p_m']:.3f}) | complex-only clusters {int(only_p.sum())} "
                  f"({100 * row['complex_mass_in_them']:.1f} % of the complex), water-only {int(only_w.sum())} "
                  f"({100 * row['water_mass_in_them']:.1f} % of water)", flush=True)
        pd.DataFrame(rows).to_csv(f"{args.out}/kl_joint_cut_{name}.csv", index=False)
        print(f"wrote {args.out}/kl_joint_cut_{name}.csv", flush=True)


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    files = sorted(f for f in glob.glob(f"{args.out}/kl_joint_cut_*.csv"))
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    if len(set(df.settings)) != 1:
        sys.exit(f"settings differ between compounds: {set(df.settings)}")
    info = ligand_info()
    df["rot_bonds"] = [info[l][0] for l in df.ligand]
    df.drop(columns="settings").to_csv(f"{args.out}/kl_joint_cut.csv", index=False)
    lo, hi = int(df.rot_bonds.min()), int(df.rot_bonds.max())
    cmap = plt.get_cmap("viridis", hi - lo + 1)
    norm = BoundaryNorm(np.arange(lo - 0.5, hi + 1.5), cmap.N)
    sub = (f"joint barrier-cut clustering of water + complex frames (equal weight), REST2 state 0; "
           f"bars: block SE over {args.blocks} time blocks")

    def colourbar(fig, ax):
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6, pad=0.02,
                          ticks=range(lo, hi + 1))
        cb.set_label("rotatable bonds")

    def scatter(ax, xcol, ycol, top, star=None):
        ax.plot([0, top], [0, top], color="0.6", ls="--", lw=1, zorder=1)
        for _, r in df.iterrows():
            ax.errorbar(r[xcol], r[ycol], xerr=r[f"{xcol}_se"], yerr=r[f"{ycol}_se"], fmt="o", ms=10,
                        mfc=cmap(norm(r.rot_bonds)), mec="k", mew=1, ecolor="0.25", elinewidth=1, capsize=3, zorder=4)
            mark = "*" if star is not None and star(r) else ""
            ax.annotate(point_label(df, r) + mark, (r[xcol], r[ycol]), textcoords="offset points", xytext=(7, 4),
                        fontsize=8.5, zorder=5)
        ax.set_xlim(-0.03 * top, top)
        ax.set_ylim(0, top)
        ax.grid(alpha=0.3)

    # 1. joint KL, both directions
    top = float(np.nanmax(np.r_[df.KL_p_w + df.KL_p_w_se, df.KL_w_p + df.KL_w_p_se])) * 1.12
    fig, ax = plt.subplots(figsize=(10.5, 8.5), layout="constrained")
    starred = lambda r: r.complex_mass_in_them > 0.01 or r.water_mass_in_them > 0.01
    scatter(ax, "KL_w_p", "KL_p_w", top, star=starred)
    ax.set_xlabel(r"KL($\pi^{water}\,\|\,\pi^{protein}$)  (k$_B$)")
    ax.set_ylabel(r"KL($\pi^{protein}\,\|\,\pi^{water}$)  (k$_B$)")
    ax.set_title(f"Bound vs free ligand ensemble on ONE joint clustering\n{sub}\n* > 1 % of one ensemble in clusters "
                 f"the other never visits: those clusters carry a 0.5-frame pseudocount", fontsize=10)
    colourbar(fig, ax)
    fig.savefig(f"{args.out}/kl_joint_cut.png", dpi=150)
    plt.close(fig)

    # 2. JS halves: KL(w||m) vs KL(p||m), both in [0, ln 2]
    fig, ax = plt.subplots(figsize=(10.5, 8.5), layout="constrained")
    scatter(ax, "JS_w_m", "JS_p_m", LN2 * 1.06)
    for f in (ax.axhline, ax.axvline):
        f(LN2, color="0.3", ls=":", lw=1.2)
    ax.text(0.01, LN2, " ln 2: no cluster shared", va="bottom", fontsize=8, color="0.3")
    ax.set_xlabel(r"KL($\pi^{water}\,\|\,m$)  (k$_B$),  m = ($\pi^{water}$ + $\pi^{protein}$)/2")
    ax.set_ylabel(r"KL($\pi^{protein}\,\|\,m$)  (k$_B$)")
    ax.set_title("Jensen-Shannon divergence in its two halves, JS = (x + y)/2\n" + sub, fontsize=10)
    colourbar(fig, ax)
    fig.savefig(f"{args.out}/js_joint_halves_cut.png", dpi=150)
    plt.close(fig)

    # 3. JS ranked
    d = df.sort_values("JS").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.5, 1.5 + 0.38 * len(d)), layout="constrained")
    ax.barh(range(len(d)), d.JS, xerr=d.JS_se, color=[cmap(norm(v)) for v in d.rot_bonds], edgecolor="k",
            lw=0.8, capsize=3, error_kw={"ecolor": "0.25", "elinewidth": 1})
    ax.set_yticks(range(len(d)), [point_label(df, r) for r in d.itertuples()])
    ax.axvline(LN2, color="0.3", ls=":", lw=1.2)
    ax.text(LN2, len(d) - 0.4, " ln 2", va="center", fontsize=8, color="0.3")
    ax.set_xlim(0, LN2 * 1.08)
    ax.set_xlabel(r"JS($\pi^{water}$, $\pi^{protein}$)  (k$_B$):  0 = same ensemble, ln 2 = no cluster shared")
    ax.grid(alpha=0.3, axis="x")
    ax.set_title("Jensen-Shannon divergence, bound vs free ligand ensemble\n" + sub, fontsize=9.5)
    colourbar(fig, ax)
    fig.savefig(f"{args.out}/js_joint_ranked_cut.png", dpi=150)
    plt.close(fig)
    print(f"wrote {args.out}/kl_joint_cut.png, js_joint_halves_cut.png, js_joint_ranked_cut.png, kl_joint_cut.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--compounds", nargs="+", default=DEFAULT)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    plot(args) if args.plot else run(args)


if __name__ == "__main__":
    main()
