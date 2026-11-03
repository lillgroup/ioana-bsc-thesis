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
  python -I kl_cross.py --plot     -> <out>/kl_cross_{cut,hdbscan}.png + .csv (all compounds)
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


def plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    titles = {"cut": "torsion partitioning (barrier cuts)", "hdbscan": "symmetry-first HDBSCAN (no single cluster)"}
    for method, title in titles.items():
        import glob
        files = sorted(glob.glob(f"{args.out}/kl_cross_{method}_*.csv"))
        if not files:
            print(f"skip {method}: no kl_cross_{method}_*.csv")
            continue
        df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        if len(set(df.settings)) != 1:
            sys.exit(f"{method}: settings differ between compounds: {set(df.settings)}")
        df.drop(columns="settings").to_csv(f"{args.out}/kl_cross_{method}.csv", index=False)
        ligands = list(dict.fromkeys(df.ligand))
        palette = list(plt.get_cmap("tab20").colors[::2]) + list(plt.get_cmap("tab20").colors[1::2]) + list(plt.get_cmap("Dark2").colors)
        colours = dict(zip(ligands, palette))
        fig, ax = plt.subplots(figsize=(10.5, 8.5), layout="constrained")
        top = np.nanmax(np.r_[df.KL_p_w_waterdef + df.KL_p_w_se, (df.KL_w_p_complexdef + df.KL_w_p_se).dropna()]) * 1.12
        ax.plot([0, top], [0, top], color="0.6", ls="--", lw=1, label="KL(w||p) = KL(p||w)")
        missing = []
        for _, r in df.iterrows():
            col = colours[r.ligand]
            if pd.isna(r.KL_w_p_complexdef):
                missing.append(f"{r.ligand} {r.complex}: KL(p||w) = {r.KL_p_w_waterdef:.2f}, no cluster in the complex")
                ax.errorbar(0, r.KL_p_w_waterdef, yerr=r.KL_p_w_se, fmt="s", mfc="white", mec=col, ecolor=col, ms=9,
                            mew=2, capsize=3, clip_on=False)
                ax.annotate(f"{r.ligand} {r.complex}\n(no complex cluster)", (0, r.KL_p_w_waterdef),
                            textcoords="offset points", xytext=(8, -4), fontsize=7.5, color=col)
                continue
            ax.errorbar(r.KL_w_p_complexdef, r.KL_p_w_waterdef, xerr=r.KL_w_p_se, yerr=r.KL_p_w_se, fmt="o",
                        color=col, ms=9, capsize=3)
            ax.annotate(f"{r.ligand} {r.complex}", (r.KL_w_p_complexdef, r.KL_p_w_waterdef),
                        textcoords="offset points", xytext=(7, 4), fontsize=8.5, color=col)
        ax.set_xlim(-0.03 * top, top)
        ax.set_ylim(0, top)
        ax.set_xlabel(r"KL($\pi^{water}\,\|\,\pi^{protein}$)  (k$_B$), clusters defined on the complex")
        ax.set_ylabel(r"KL($\pi^{protein}\,\|\,\pi^{water}$)  (k$_B$), clusters defined in water")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=8.5)
        ax.set_title(f"Bound vs free ligand ensemble, clusters of one ensemble, frames of the other assigned\n"
                     f"{title}, REST2 state 0; bars: block SE over {args.blocks} time blocks", fontsize=10)
        fig.savefig(f"{args.out}/kl_cross_{method}.png", dpi=150)
        plt.close(fig)
        print(f"wrote {args.out}/kl_cross_{method}.png" + (f" ({'; '.join(missing)})" if missing else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--method", choices=("cut", "hdbscan"))
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--compounds", nargs="+", default=DEFAULT)
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    if args.plot:
        plot(args)
    elif args.method:
        run(args)
    else:
        ap.error("--method or --plot")


if __name__ == "__main__":
    main()
