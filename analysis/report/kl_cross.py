"""KL divergence between the bound and the free ligand ensemble with Shu-Yu's protocol: clusters
defined on ONE ensemble, the other ensemble's frames assigned to them, in both directions.

  water-defined clusters:   cluster the water ensemble, assign the complex frames
                            -> pi^w, pi^p on the water clusters ->  KL(pi^p || pi^w)   (y axis)
  complex-defined clusters: cluster the complex ensemble, assign the water frames
                            -> pi^p, pi^w on the complex clusters -> KL(pi^w || pi^p)  (x axis)

The reference ensemble defined the clusters, so its populations are > 0 everywhere and both KL are
finite without a pseudocount (sum_c p ln(p/q), p = assigned ensemble, q = reference).

  --method cut      torsion partitioning (cut_clustering, MD-tools 306b109, env analysis-env-cut),
                    defaults, every frame. MD-tools has no assignment for it: a frame is walked down
                    the cut tree with the rule _Cutter uses to split (periodic cut c1, c2: first child
                    iff (x - c1) mod 360 < (c2 - c1) mod 360; bounded cut c in a cell from lo: iff
                    (x - lo) mod 360 < (c - lo) mod 360), on the symmetry-invariant coordinates computed
                    by MD-tools' own _invariant_coordinates (names checked against the reference).
  --method hdbscan  symmetry-first t_hdbscan, allow_single_cluster=False (MD-tools c30c1c0, env
                    analysis-env-main), fitted on the reference at --stride (<= 10 000 frames), the
                    other ensemble assigned by THDBSCAN.predict at the same stride. Noise / ambiguous
                    frames are left out and populations renormalised; the assigned ensemble's
                    unassigned fraction is reported (complex frames no water cluster takes: induced fit).
                    A reference with no cluster (NO_CLUSTER, e.g. 4ibi bound) gives no value.

Both ensembles are taken in the water ligand's atom order (the complex reordered by prepare.atom_map;
folate 3QL3's tautomeric H takes the water's H slot, which no torsion uses), torsions from the water
ligand's files (prepare.py, one numbering). Error bars: block SE over --blocks time blocks of the assigned
ensemble against the reference's full populations (the clustering is not refitted).

  python -I kl_cross.py --method cut|hdbscan --compounds 4ibi  -> <out>/kl_cross_<method>_4ibi.csv
  (kl_cross.sbatch: one Slurm task per compound and method, then --plot)
  python -I kl_cross.py --plot     -> <out>/kl_cross_{cut,hdbscan}.png + .csv (all compounds), and
                                      <out>/kl_cross_{cut,hdbscan}_regions.png: the same with the binding
                                      regimes as background (regime(); --lock-key-max, --induced-fit-min);
                                      <out>/kl_cross_cut_examples{,_regions}.png: only 4j7d, 4crf, 4ase, 1l8g.
                                      Points coloured by rotatable bonds
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
import csv as _csv
DEFAULT = [r["name"] for r in _csv.DictReader(open(os.path.join(HERE, "compounds.csv")))
           if r["complexes"]]                                          # paracetamol has no complex
KINDS = ("torsions", "torsions_rigid", "torsions_polarH")


def kl(p, q):
    m = p > 0
    return float((p[m] * np.log(p[m] / q[m])).sum())


def pops(labels, k):
    lab = labels[(labels >= 0) & (labels < k)]
    return np.bincount(lab, minlength=k) / max(len(lab), 1)


def ensembles(man, c):
    """water and complex trajectories, both in the water ligand's atom order; water quads and names."""
    import mdtraj as md
    from rdkit import Chem
    from cluster_prune import load_torsions
    from prepare import atom_map
    w = man["water"]
    tw = md.load(w["traj"], top=w["top"])
    tc = md.load(c["traj"], top=c["top"])
    amap, _ = atom_map(Chem.MolFromMolFile(w["sdf"], removeHs=False), Chem.MolFromMolFile(c["sdf"], removeHs=False))
    spare = iter(sorted(set(range(tc.n_atoms)) - {v for v in amap.values() if v is not None}))
    order = [amap[i] if amap.get(i) is not None else next(spare) for i in range(tw.n_atoms)]
    tc = md.Trajectory(tc.xyz[:, order], tw.topology, time=tc.time)
    tors = load_torsions([f"{w['torsions']}/{k}.txt" for k in KINDS if os.path.getsize(f"{w['torsions']}/{k}.txt")])
    return tw, tc, list(tors), list(tors.values())


def thdbscan(traj, mol, names, quads, extra, **kw):
    import mdtraj as md
    from md_tools.analysis import t_hdbscan, TorsionDefinitions, verify_atom_mapping
    mapping, _ = verify_atom_mapping(mol, [(a.index, b.index) for a, b in traj.topology.bonds],
                                     [a.element.symbol for a in traj.topology.atoms])
    defs = TorsionDefinitions(quads=tuple(quads), names=tuple(names), indexing="zero-based", units="radian")
    return t_hdbscan(md.compute_dihedrals(traj, quads), molecule=mol, torsion_definitions=defs,
                     coordinates=traj, atom_mapping=mapping, **extra, **kw)


def cut_assign(r, Y):
    """Leaf cluster of every row of Y (frames x cut coordinates) in the cut tree of result r."""
    tree = {n["node"]: n for n in r["tree"]}
    col = {name: j for j, name in enumerate(r["coordinates"])}
    root = next(i for i, n in tree.items() if n["parent"] is None)
    out = np.empty(len(Y), dtype=np.int64)
    for f, y in enumerate(np.degrees(Y)):
        i = root
        while tree[i]["children"]:
            cut = tree[i]["cut"]
            x = y[col[cut["coordinate"]]]
            if cut["kind"] == "periodic":
                c1, c2 = cut["positions_degrees"]
                first = (x - c1) % 360.0 < (c2 - c1) % 360.0
            else:
                lo = tree[i]["box_degrees"][cut["coordinate"]][0]
                first = (x - lo) % 360.0 < (cut["positions_degrees"][0] - lo) % 360.0
            i = tree[i]["children"][0 if first else 1]
        out[f] = tree[i]["cluster"]
    return out


class Ensemble:
    """One state-0 ensemble and what the methods need from it, each built ONCE (MD-tools enumerates the
    molecule's symmetry whenever a t_hdbscan object is made -- minutes for 4j7d -- so the objects are
    reused in both directions and, for folate, for both complexes)."""

    def __init__(self, label, traj, method, mol, names, quads, extra, stride):
        self.label, self.traj, self.method = label, traj, method
        self.mol, self.names, self.quads, self.extra, self.stride = mol, names, quads, extra, stride
        self._cut = self._fit = self._Y = None
        self.no_cluster = False

    def cut(self):                      # cut_clustering of this ensemble (as reference)
        if self._cut is None:
            from md_tools.analysis import cut_clustering
            from md_tools.analysis._t_cut import _invariant_coordinates
            obj = thdbscan(self.traj, self.mol, self.names, self.quads, self.extra)
            self._cut = cut_clustering(obj)
            self._Y = _invariant_coordinates(obj, self._cut["settings"]["z_min"])[:2]
        return self._cut

    def coordinates(self):              # MD-tools' cut coordinates of this ensemble (as query)
        self.cut()
        return self._Y

    def fit(self):                      # HDBSCAN fit of this ensemble (as reference), or None
        if self._fit is None and not self.no_cluster:
            try:
                self._fit = thdbscan(self.traj[::self.stride], self.mol, self.names, self.quads, self.extra,
                                     allow_single_cluster=False).fit()
            except RuntimeError as e:   # one-state ensemble (cluster_prune's NO_CLUSTER case)
                if "no draw produced a cluster" not in str(e):
                    raise
                self.no_cluster = True
        return self._fit


def one_direction(ref, qry, blocks):
    """Clusters of `ref`, frames of `qry` assigned: KL(pi_qry || pi_ref), block SE, ... or None."""
    if ref.method == "cut":
        r = ref.cut()
        k, lr = r["n_clusters"], r["labels"]
        Yq, cq = qry.coordinates()
        if list(cq) != list(r["coordinates"]):
            sys.exit(f"cut coordinates differ: {cq} vs {r['coordinates']}")
        lq = cut_assign(r, Yq)
        settings = {k_: v for k_, v in r["settings"].items() if k_ not in ("floor_frames", "reference")}
    else:
        import mdtraj as md
        fit = ref.fit()
        if fit is None:
            return None
        k, lr = fit.n_clusters_, fit.labels_
        q = qry.traj[::qry.stride]
        lq, _ = fit.predict(md.compute_dihedrals(q, qry.quads), coordinates=q)
        s = fit.summary()
        settings = {"route": s["route"], "allow_single_cluster": s["allow_single_cluster"], "mass_floor": s["mass_floor"]}
    pr, pq = pops(lr, k), pops(lq, k)
    # block SE: the assigned ensemble in time blocks against the reference's populations over ALL its
    # frames (the reference defines the clusters; per-block reference populations can be 0 for a small
    # cluster and make one block's KL explode)
    vals = [kl(pops(bq, k), pr) for bq in np.array_split(lq, blocks)]
    return {"KL": kl(pq, pr), "se": float(np.std(vals, ddof=1) / np.sqrt(len(vals))), "k": int(k),
            "unassigned": float(np.mean((lq < 0) | (lq >= k))), "ref_pops": pr, "qry_pops": pq,
            "settings": settings}


def run(args):
    import md_tools
    from rdkit import Chem
    os.makedirs(args.out, exist_ok=True)
    for name in args.compounds:
        rows = []
        man = json.load(open(f"{WORK}/{name}/manifest.json"))
        mol = Chem.SDMolSupplier(man["water"]["sdf"], removeHs=False)[0]
        extra = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
        stride = max(man["stride"], 1)
        water = None
        for c in man["complexes"]:
            tw, tc, names, quads = ensembles(man, c)
            if water is None:
                water = Ensemble("water", tw, args.method, mol, names, quads, extra, stride)
            cplx = Ensemble(c["label"], tc, args.method, mol, names, quads, extra, stride)
            wdef = one_direction(water, cplx, args.blocks)       # water clusters, complex frames assigned
            pdef = one_direction(cplx, water, args.blocks)       # complex clusters, water frames assigned
            if wdef is None:
                sys.exit(f"{name}: no cluster in water")
            row = {"method": args.method, "ligand": name, "complex": c["label"],
                   "KL_p_w_waterdef": round(wdef["KL"], 4), "KL_p_w_se": round(wdef["se"], 4),
                   "n_water_clusters": wdef["k"], "complex_unassigned_on_water_clusters": round(wdef["unassigned"], 4),
                   "top_complex_share_on_water_clusters": round(float(wdef["qry_pops"].max()), 4),
                   "water_share_of_that_cluster": round(float(wdef["ref_pops"][wdef["qry_pops"].argmax()]), 4),
                   "KL_w_p_complexdef": None if pdef is None else round(pdef["KL"], 4),
                   "KL_w_p_se": None if pdef is None else round(pdef["se"], 4),
                   "n_complex_clusters": 0 if pdef is None else pdef["k"],
                   "water_unassigned_on_complex_clusters": None if pdef is None else round(pdef["unassigned"], 4),
                   "note": "" if pdef is not None else "no cluster in the complex: complex-defined KL undefined",
                   "md_tools": md_tools.__version__, "settings": json.dumps(wdef["settings"], sort_keys=True, default=str)}
            rows.append(row)
            print(f"{name:11s} {c['label']}: water-defined {wdef['k']:3d} clusters  KL(p||w) {wdef['KL']:.2f} +- {wdef['se']:.2f}"
                  f" (complex unassigned {100 * wdef['unassigned']:.1f} %) | complex-defined "
                  + ("-- no cluster" if pdef is None else
                     f"{pdef['k']:3d} clusters  KL(w||p) {pdef['KL']:.2f} +- {pdef['se']:.2f} (water unassigned "
                     f"{100 * pdef['unassigned']:.1f} %)"), flush=True)
        pd.DataFrame(rows).to_csv(f"{args.out}/kl_cross_{args.method}_{name}.csv", index=False)
        print(f"wrote {args.out}/kl_cross_{args.method}_{name}.csv", flush=True)


EXAMPLES = ("4j7d", "4crf", "4ase", "1l8g")     # the four corners of the cut plot (Shu-Yu: explain the equation with them)
LABEL_OFFSETS = {"5ov8": (-9, 5, "right"), "4p6w": (8, -12, "left"), "1eld": (8, 5, "left")}  # crowded corner
SETS = {"folate": "^", "original small molecules": "s", "rotatable-bond set": "o"}


def ligand_info():
    """name -> (rotatable bonds, compound set): rotatable bonds of the water ligand as in
    conformational_entropy.py (RDKit, heavy atoms), set from compounds.csv."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    out = {}
    for row in _csv.DictReader(open(os.path.join(HERE, "compounds.csv"))):
        path = f"{WORK}/{row['name']}/manifest.json"
        if not os.path.exists(path):
            continue
        sdf = json.load(open(path))["water"]["sdf"]
        rb = rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(Chem.MolFromMolFile(sdf)))
        kind = ("folate" if row["name"] == "folic-acid" else
                "rotatable-bond set" if row["water"].startswith("REST2-rotbonds") else "original small molecules")
        out[row["name"]] = (int(rb), kind)
    return out


def point_label(df, r):
    """The ligand name; the complex only where a ligand has several (folate 4LRH / 3QL3)."""
    return f"{r.ligand} {r.complex}" if (df.ligand == r.ligand).sum() > 1 else r.ligand


def draw_kl(ax, df, info, cmap, norm, top, args, regions, big=False):
    from matplotlib.lines import Line2D
    if regions:
        draw_regions(ax, top, args)
    ax.plot([0, top], [0, top], color="0.6", ls="--", lw=1, zorder=1)
    for _, r in df.iterrows():
        rb, kind = info[r.ligand]
        col = cmap(norm(rb))
        x = 0 if pd.isna(r.KL_w_p_complexdef) else r.KL_w_p_complexdef
        ax.errorbar(x, r.KL_p_w_waterdef, xerr=None if pd.isna(r.KL_w_p_complexdef) else r.KL_w_p_se,
                    yerr=r.KL_p_w_se, fmt="o", ms=15 if big else 10, mfc="white" if pd.isna(r.KL_w_p_complexdef) else col,
                    mec="k" if not pd.isna(r.KL_w_p_complexdef) else col, mew=1.0 if not pd.isna(r.KL_w_p_complexdef) else 2,
                    ecolor="0.25", elinewidth=1, capsize=3, zorder=4, clip_on=False)
        text = point_label(df, r) + ("\n(no complex cluster)" if pd.isna(r.KL_w_p_complexdef) else "")
        dx, dy, ha = (10, 6, "left") if big else LABEL_OFFSETS.get(r.ligand, (7, 4, "left"))
        ax.annotate(text, (x, r.KL_p_w_waterdef), textcoords="offset points", xytext=(dx, dy), ha=ha,
                    fontsize=13 if big else 8.5, color="k", fontweight="bold" if big else "normal", zorder=5)
    ax.set_xlim(-0.03 * top, top)
    ax.set_ylim(0, top)
    ax.set_xlabel(r"KL($\pi^{water}\,\|\,\pi^{protein}$)  (k$_B$), clusters defined on the complex")
    ax.set_ylabel(r"KL($\pi^{protein}\,\|\,\pi^{water}$)  (k$_B$), clusters defined in water")
    ax.grid(alpha=0.3)
    handles = [Line2D([], [], color="0.6", ls="--", lw=1, label="KL(w||p) = KL(p||w)")]
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, framealpha=0.9)


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    titles = {"cut": "torsion partitioning (barrier cuts)", "hdbscan": "symmetry-first HDBSCAN (no single cluster)"}
    info = ligand_info()
    for method, title in titles.items():
        import glob
        files = sorted(glob.glob(f"{args.out}/kl_cross_{method}_*.csv"))
        files = [f for f in files if not f.endswith(("_regions.csv", "_examples.csv"))]
        if not files:
            print(f"skip {method}: no kl_cross_{method}_*.csv")
            continue
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        if len(set(df.settings)) != 1:
            sys.exit(f"{method}: settings differ between compounds: {set(df.settings)}")
        df["regime"] = [regime(x, y, args) for x, y in zip(df.KL_w_p_complexdef, df.KL_p_w_waterdef)]
        df["rot_bonds"] = [info[l][0] for l in df.ligand]
        df["compound_set"] = [info[l][1] for l in df.ligand]
        df.drop(columns="settings").to_csv(f"{args.out}/kl_cross_{method}.csv", index=False)
        top = np.nanmax(np.r_[df.KL_p_w_waterdef + df.KL_p_w_se, (df.KL_w_p_complexdef + df.KL_w_p_se).dropna()]) * 1.12
        lo, hi = int(df.rot_bonds.min()), int(df.rot_bonds.max())
        cmap = plt.get_cmap("viridis", hi - lo + 1)
        norm = BoundaryNorm(np.arange(lo - 0.5, hi + 1.5), cmap.N)
        sub = f"{title}, REST2 state 0; bars: block SE over {args.blocks} time blocks"
        variants = [("", df, False), ("_regions", df, True)]
        if method == "cut":
            ex = df[df.ligand.isin(EXAMPLES)]
            variants += [("_examples", ex, False), ("_examples_regions", ex, True)]
        for suffix, d, regions in variants:
            big = "examples" in suffix
            fig, ax = plt.subplots(figsize=(10.5 + 1.2 * regions, 8.5), layout="constrained")
            draw_kl(ax, d, info, cmap, norm, top, args, regions, big=big)
            cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6,
                              pad=0.09 if regions else 0.02, ticks=range(lo, hi + 1))
            cb.set_label("rotatable bonds")
            ax.set_title(f"Bound vs free ligand ensemble, clusters of one ensemble, frames of the other assigned\n{sub}"
                         + ("\nunassigned (noise) frames are dropped before the KL, so regimes are unreliable where "
                            "most frames are unassigned: use the barrier-cut plot" if regions and method == "hdbscan" else ""),
                         fontsize=10)
            path = f"{args.out}/kl_cross_{method}{suffix}.png"
            fig.savefig(path, dpi=150)
            plt.close(fig)
            print(f"wrote {path}")
        print("  " + "; ".join(f"{g}: {', '.join(point_label(df, r) for r in d.itertuples())}"
                               for g, d in df.groupby("regime", sort=False)))


REGIMES = {"lock and key": "#2ca25f", "conformational selection": "#3182bd", "induced fit": "#e6550d"}


def regime(x, y, args):
    """Binding regime from the two KL (k_B): y = KL(p||w) ~ -ln(population of the bound conformer in
    water) when the complex sits in one cluster; x = KL(w||p) = how differently the free ligand fills
    the complex's clusters.
      induced fit               y >= --induced-fit-min: the bound conformer is (nearly) absent in water
      lock and key              x and y < --lock-key-max: free and bound ensembles nearly the same
      conformational selection  otherwise: the bound conformer is a populated state of the free
                                ligand, the protein shifts the populations towards it
    With no complex cluster (x undefined) a y below --lock-key-max cannot be classified."""
    if y >= args.induced_fit_min:
        return "induced fit"
    if pd.isna(x):
        return "conformational selection" if y >= args.lock_key_max else "unclassified"
    if max(x, y) < args.lock_key_max:
        return "lock and key"
    return "conformational selection"


def draw_regions(ax, top, args):
    """Background of the three regimes (regime()) and a right axis e^-y: the bound conformer's
    approximate population in water."""
    from matplotlib.patches import Patch, Rectangle
    lk, fit = args.lock_key_max, args.induced_fit_min
    ax.add_patch(Rectangle((-1, fit), top + 2, top, color=REGIMES["induced fit"], alpha=0.13, lw=0, zorder=0))
    ax.add_patch(Rectangle((-1, 0), top + 2, fit, color=REGIMES["conformational selection"], alpha=0.13, lw=0, zorder=0))
    ax.add_patch(Rectangle((-1, 0), lk + 1, lk, color="white", lw=0, zorder=0))
    ax.add_patch(Rectangle((-1, 0), lk + 1, lk, color=REGIMES["lock and key"], alpha=0.22, lw=0, zorder=0))
    rules = {"lock and key": f"KL(w||p), KL(p||w) < {lk:g}: free ≈ bound ensemble",
             "conformational selection": f"bound conformer populated in water (KL(p||w) < {fit:g}, e$^{{-{fit:g}}}$ ≈ "
                                         f"{100 * np.exp(-fit):.0f} %), protein shifts populations",
             "induced fit": f"KL(p||w) ≥ {fit:g}: bound conformer ≲ {100 * np.exp(-fit):.0f} % of the free ensemble"}
    for name, (x, y) in {"lock and key": (0.52 * lk, 0.01 * top), "conformational selection": (0.97 * top, 0.04 * top),
                         "induced fit": (0.97 * top, fit + 0.04 * top)}.items():
        ax.text(x, y, name, color=REGIMES[name], fontsize=11 if name == "lock and key" else 13, fontweight="bold",
                alpha=0.75,
                ha="left" if name == "lock and key" else "right", va="bottom", zorder=1)
    leg = ax.legend(handles=[Patch(color=REGIMES[n], alpha=0.35, label=f"{n}: {rules[n]}") for n in REGIMES],
                    loc="upper right", bbox_to_anchor=(0.99, 0.82), fontsize=8, title="rules", title_fontsize=8.5,
                    framealpha=0.9)
    ax.add_artist(leg)
    right = ax.twinx()
    pcts = [p for p in (50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01) if -np.log(p / 100) <= top]
    right.set_ylim(0, top)
    right.set_yticks([-np.log(p / 100) for p in pcts], [f"{p:g} %" for p in pcts])
    right.set_ylabel("≈ population of the bound conformer in water, e$^{-KL(p||w)}$\n"
                     "(when the complex sits mostly in one cluster)", fontsize=9)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--method", choices=("cut", "hdbscan"))
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--compounds", nargs="+", default=DEFAULT)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--lock-key-max", type=float, default=1.0,
                    help="regions plot: lock and key below this KL (k_B) on both axes")
    ap.add_argument("--induced-fit-min", type=float, default=3.0,
                    help="regions plot: induced fit from this KL(p||w) (k_B) up (e^-3 = 5 %% in water)")
    args = ap.parse_args()
    if args.plot:
        plot(args)
    elif args.method:
        run(args)
    else:
        ap.error("--method or --plot")


if __name__ == "__main__":
    main()
