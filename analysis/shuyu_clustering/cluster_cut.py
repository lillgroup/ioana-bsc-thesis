#!/usr/bin/env python
"""
Shu-Yu's TORSION PARTITIONING (torsional clustering by barrier cuts, MD-tools main 306b109,
`cut_clustering`), with the figures of this repository.

The method, exactly as MD-tools implements it (docs: openmm_methods/analysis/torsional-clustering/
torsion-partitioning, torsion-preprocessing; tutorials ALA/analysis/clustering and
paracetamol/analysis/clustering):
  1. Coordinates   t_hdbscan(tors, molecule, torsion_definitions, coordinates, atom_mapping), NOT
                   fitted: only its angles and symmetry operations are read. A torsion no symmetry
                   operation moves is its own coordinate; torsions an operation moves together are
                   replaced by symmetry-invariant averages (e.g. 2*aryl, phenol-aryl).
  2. Empty-bin cuts, recursively, wherever both sides keep >= mass_floor of ALL frames.
  3. Annealing     each cell from E* = its largest finite free-energy range down to kT_cut in
                   kT_step steps; at each E* it is cut on the coordinate with the highest barrier
                   that clears E* (F = -ln p on n_bins bins; cuts at bin midpoints).
  4. Numbering     cells that are not cut again are the clusters, by mass. EVERY frame is labelled:
                   no noise, no ambiguous code.
Defaults are MD-tools' (kT_cut 1, mass_floor 0.001, min_p2p pi/3, kT_step 0.5, 72 bins, z_min 0.5);
it needs no distance matrix, so it runs on every frame (--stride 1 by default).

Figures and files (--out DIR):
  clusters_overview.png   ../report-style overview (cluster_prune.draw_overview; gold frame: the cluster
                          the crystal falls in, by walking its cut coordinates down the tree): molecule, one column
                          per torsion t0.., the full distribution and one row per cluster (clusters
                          >= --draw-min-mass; the rest summed in the title), all in ONE labelling
                          (result["aligned_angles"], the labelling nearest frame `reference`), crystal
                          lines in that labelling
  split_tree.png          MD-tools plot_split_tree: every cell, its 1-D profiles, the cuts
  clusters_free_energy.png MD-tools plot_hist(kind="free_energy") for the largest clusters
  cut_log.txt, groups.csv, labels.csv, summary.json
Torsion numbering: the order of the --torsions-file files and lines, as in cluster_prune.py, so t<k>
is the same bond in every figure of a compound (../report/prepare.py writes one set per system).

  python cluster_cut.py --traj ... --top ... --sdf ... --torsions-file torsions.txt \
      [--torsions-file torsions_rigid.txt ...] [--ref-sdf crystal.sdf] --name "4ibi water" --out DIR
Needs the environment with MD-tools >= 306b109 (analysis-env-cut).
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_prune import (load_torsions, draw_overview, circular_mean, circular_std,  # noqa: E402
                           read_rest2_config, simulated_ns)

try:
    from md_tools.analysis import (t_hdbscan, cut_clustering, plot_split_tree, plot_hist,
                                   TorsionDefinitions, verify_atom_mapping, group_name,
                                   reference_assignments)
except ImportError as e:
    sys.exit(f"md_tools.analysis has no cut_clustering ({e}): needs MD-tools main >= 306b109 "
             "(env analysis-env-cut)")


def cut_assign(r, Y):
    """Leaf cluster of every row of Y (rows x cut coordinates, radians) in the cut tree of result r,
    or -1 where a cut is made on a coordinate the row does not define (NaN). The side of a cut is
    decided as MD-tools' _Cutter splits: periodic cut c1, c2 -> first child iff (x - c1) mod 360 <
    (c2 - c1) mod 360; bounded cut c in a cell from lo -> iff (x - lo) mod 360 < (c - lo) mod 360.
    (A copy lives in ../report/kl_cross.py, which must also run in the env without cut_clustering.)"""
    tree = {nd["node"]: nd for nd in r["tree"]}
    col = {name: j for j, name in enumerate(r["coordinates"])}
    root = next(i for i, nd in tree.items() if nd["parent"] is None)
    out = np.empty(len(Y), dtype=np.int64)
    for f, y in enumerate(np.degrees(Y)):
        i = root
        while tree[i]["children"]:
            cut = tree[i]["cut"]
            x = y[col[cut["coordinate"]]]
            if not np.isfinite(x):
                i = None
                break
            if cut["kind"] == "periodic":
                c1, c2 = cut["positions_degrees"]
                first = (x - c1) % 360.0 < (c2 - c1) % 360.0
            else:
                lo = tree[i]["box_degrees"][cut["coordinate"]][0]
                first = (x - lo) % 360.0 < (cut["positions_degrees"][0] - lo) % 360.0
            i = tree[i]["children"][0 if first else 1]
        out[f] = -1 if i is None else tree[i]["cluster"]
    return out


def crystal_cluster(fit, r, rows):
    """The cluster a crystal ligand falls in: its cut coordinates (MD-tools' _invariant_coordinates on
    the crystal's angle rows -- every atom assignment, which the invariants do not distinguish) walked
    down the cut tree. Returns (cluster or None, note)."""
    from types import SimpleNamespace
    from md_tools.analysis._t_cut import _invariant_coordinates
    proxy = SimpleNamespace(n_torsions_=fit.n_torsions_, names=fit.names, symmetry=fit.symmetry,
                            descriptor_=fit.descriptor_, closed_theta_=np.asarray(rows, float))
    with np.errstate(invalid="ignore"):
        Y, cname, _, _ = _invariant_coordinates(proxy, r["settings"]["z_min"])
    if list(cname) != list(r["coordinates"]):
        return None, f"cut coordinates differ ({cname})"
    found = set(int(c) for c in cut_assign(r, Y))
    if found == {-1} or -1 in found and len(found) == 1:
        return None, "a cut on its way is on a coordinate the crystal does not define (through H)"
    found.discard(-1)
    if len(found) > 1:
        return None, f"its atom assignments fall in different clusters ({', '.join(group_name(c) for c in sorted(found))})"
    return found.pop(), ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--traj", required=True, help="state-0 trajectory of the ligand (or solute)")
    ap.add_argument("--top", required=True, help="matching topology WITH bonds")
    ap.add_argument("--sdf", required=True, help="ligand SDF (bond orders), the trajectory's atom order")
    ap.add_argument("--torsions-file", action="append", required=True)
    ap.add_argument("--stride", type=int, default=1, help="default 1: the method needs no distance matrix")
    ap.add_argument("--ref-sdf", action="append", default=None, help="crystal ligand(s): lines in the figures")
    ap.add_argument("--name", default=None)
    ap.add_argument("--rest2-config", default=None)
    ap.add_argument("--max-matches", type=int, default=None, help="cap on graph automorphisms (as cluster_prune)")
    ap.add_argument("--kT-cut", type=float, default=None, help="lowest barrier cut, kT (MD-tools default 1)")
    ap.add_argument("--mass-floor", type=float, default=None, help="smallest piece, fraction of ALL frames (default 0.001)")
    ap.add_argument("--min-p2p", type=float, default=None, help="smallest cut-to-cut distance, rad (default pi/3)")
    ap.add_argument("--draw-min-mass", type=float, default=0.01,
                    help="overview: draw clusters with at least this mass (default 1 %%); drawing only")
    ap.add_argument("--hist-clusters", type=int, default=8, help="clusters drawn by plot_hist (largest first)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    P = lambda f: os.path.join(args.out, f)

    import md_tools
    print(f"md-tools {md_tools.__version__} from {os.path.dirname(md_tools.__file__)}")
    torsions = load_torsions(args.torsions_file)
    names, quads = list(torsions), list(torsions.values())
    cmd = md.load(args.traj, top=args.top, stride=args.stride)
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    if mol is None:
        sys.exit(f"RDKit could not read {args.sdf}")
    bonds = [(b[0].index, b[1].index) for b in cmd.topology.bonds]
    if not bonds:
        sys.exit("the topology has no bonds - use one with CONECT records / bonds")
    mapping, why = verify_atom_mapping(mol, bonds, [a.element.symbol for a in cmd.topology.atoms])
    tors = md.compute_dihedrals(cmd, quads)
    definitions = TorsionDefinitions(quads=tuple(quads), names=tuple(names), indexing="zero-based", units="radian")
    extra = {} if args.max_matches is None else {"max_matches": int(args.max_matches)}
    for k, (name, q) in enumerate(torsions.items()):
        print(f"  t{k} {name} {q}")
    print(f"{tors.shape[0]} frames (stride {args.stride}), {len(names)} torsions | atom mapping: {why}")

    # ------------------------------------------------------------------ 1. the partition
    fit = t_hdbscan(tors, molecule=mol, torsion_definitions=definitions, coordinates=cmd,
                    atom_mapping=mapping, **extra)                          # NOT fitted (Shu-Yu)
    kw = {k: v for k, v in (("kT_cut", args.kT_cut), ("mass_floor", args.mass_floor),
                            ("min_p2p", args.min_p2p)) if v is not None}
    r = cut_clustering(fit, **kw)
    labels, n = r["labels"], r["n_clusters"]
    masses = [r["cluster_mass"][c] for c in range(n)]
    print(f"\ncoordinates: {r['coordinates']} | coupling groups: {r['coupling_groups']}")
    print(f"settings: {r['settings']}")
    print("\n".join(r["log"]))
    print(f"\n{n} clusters; masses (%): {[round(100 * m, 2) for m in masses]}")
    open(P("cut_log.txt"), "w").write("\n".join(r["log"]) + "\n")
    with open(P("labels.csv"), "w") as f:
        f.write("frame,time_ps,cluster\n")
        for i, (t_ps, lab) in enumerate(zip(cmd.time, labels)):
            f.write(f"{i * args.stride},{float(t_ps):.1f},{group_name(int(lab))}\n")
    al = r["aligned_angles"]                                     # frames x torsions, one labelling
    stats = [{"cluster": group_name(c), "frames": int((labels == c).sum()), "mass": round(masses[c], 5)}
             for c in range(n)]
    json.dump({"method": "torsion partitioning (cut_clustering)", "md_tools": md_tools.__version__,
               "frames": int(len(labels)), "stride": args.stride, "n_clusters": int(n),
               "cluster_mass": {group_name(c): masses[c] for c in range(n)},
               "coordinates": r["coordinates"], "coupling_groups": r["coupling_groups"],
               "settings": r["settings"], "cuts": r["cuts"]},
              open(P("summary.json"), "w"), indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))

    # ------------------------------------------------------------------ 2. one labelling per cluster (drawing only)
    # cut_clustering's aligned_angles write every frame in the labelling nearest ONE frame (`reference`);
    # a cluster far from it can then show one conformer in two labellings (paracetamol B: half a turn).
    # As cluster_prune's align_all(per_cluster=True): each cluster is written in the labelling nearest
    # its own central frame -- the member closest to the cluster's circular mean in the cut coordinates,
    # which are labelling-independent -- and its crystal lines in that labelling. The cells are untouched.
    d = fit.descriptor_
    C = np.asarray(fit.closed_theta_, float)
    K = len(names)
    wrap = lambda x: np.angle(np.exp(1j * x))
    crystal_rows = {os.path.splitext(os.path.basename(p))[0]: reference_assignments(p, mol, d.quads)
                    for p in (args.ref_sdf or [])}

    def crystal_lines(ref_row):
        out = {}
        for label, rows_c in crystal_rows.items():
            near = d.nearest_labelling(rows_c, ref_row)["angles"]
            out[label] = [[float(wrap(a))] if np.isfinite(a) else [] for a in near[:K]]
        return out

    centre, cluster_al = {}, {}
    for c in range(n):
        members = np.flatnonzero(labels == c)
        Y = r["coordinate_values"][members]
        mean = np.angle(np.exp(1j * Y).mean(axis=0))
        ref = int(members[np.argmin(np.abs(wrap(Y - mean)).sum(axis=1))])
        centre[c] = ref
        cluster_al[c] = wrap(d.align_to_reference(C[members], C[ref])["angles"][:, :K])
    top_refs = crystal_lines(C[int(r["settings"]["reference"])])

    # ------------------------------------------------------------------ 3. figures
    shown = [c for c in range(n) if masses[c] >= args.draw_min_mass]
    rest = 1.0 - sum(masses[c] for c in shown)
    # crystal -> cluster (gold frame, as cluster_prune's assign_reference in the HDBSCAN figures)
    best, crystal_notes, assignment = {}, [], {}
    for label, rows_c in crystal_rows.items():
        c, why = crystal_cluster(fit, r, rows_c)
        assignment[label] = {"cluster": None if c is None else group_name(c), "note": why}
        if c is None:
            crystal_notes.append(f"crystal {label}: not placed ({why})")
        elif c in shown:
            best.setdefault(group_name(c), []).append(label)
        else:
            crystal_notes.append(f"crystal {label}: in cluster {group_name(c)} ({100 * masses[c]:.2f} %, not drawn)")
        print(f"crystal {label}: cluster {assignment[label]['cluster'] or 'none'} {why}")
    summary_path = P("summary.json")
    summary = json.load(open(summary_path))
    summary["crystal_assignment"] = assignment
    json.dump(summary, open(summary_path, "w"), indent=2)
    aligned = {group_name(c): cluster_al[c] for c in shown}
    cluster_refs = {group_name(c): crystal_lines(C[centre[c]]) for c in shown}
    for row, c in zip(stats, range(n)):         # per-cluster statistics in the cluster's own labelling
        row["central_frame"] = centre[c] * args.stride
        for k in range(K):
            row[f"t{k}_mean_deg"] = round(np.rad2deg(circular_mean(cluster_al[c][:, k])), 1)
            row[f"t{k}_std_deg"] = round(np.rad2deg(circular_std(cluster_al[c][:, k])), 1)
    pd.DataFrame(stats).to_csv(P("groups.csv"), index=False)
    n_rep, tau_max = read_rest2_config(args.rest2_config)
    info = [f"{simulated_ns(cmd, args.stride):g} ns", f"{n_rep} replicas" if n_rep else None,
            f"τ_max {tau_max}" if tau_max is not None else None, f"stride {args.stride}, {len(labels)} frames"]
    s = r["settings"]
    title = (f"{args.name or 'ligand'} — REST2 " + ", ".join(x for x in info if x) +
             f"\ntorsion partitioning (barrier cuts): {n} clusters at kT_cut {s['kT_cut']:g}, mass_floor "
             f"{s['mass_floor']:g} ({s['floor_frames']:g} frames); every frame labelled"
             f"\ndrawn: the {len(shown)} clusters >= {100 * args.draw_min_mass:g} %"
             + (f" (the other {n - len(shown)} hold {100 * rest:.1f} %)" if len(shown) < n else "")
             + ("; black lines = crystal (bioactive) conformation" if crystal_rows else "")
             + ("; gold frame = cluster the crystal falls in" if best else "")
             + ("\n" + "; ".join(crystal_notes) if crystal_notes else ""))
    draw_overview(P("clusters_overview.png"), mol=mol, names=names, quads=quads, tors=al, aligned=aligned,
                  masses=[masses[c] for c in shown], noise=0.0, refs=top_refs, cluster_refs=cluster_refs,
                  top_refs=top_refs, top_label="all frames,\none labelling", title=title, best=best)
    try:
        plot_split_tree(r)
        plt.gcf().savefig(P("split_tree.png"), dpi=130, bbox_inches="tight")
    except Exception as e:                       # a one-cell result has no tree to draw
        print(f"split tree not drawn: {e}")
    plt.close("all")
    plot_hist(r, n_cluster=min(n, args.hist_clusters), kind="free_energy")
    plt.gcf().savefig(P("clusters_free_energy.png"), dpi=130, bbox_inches="tight")
    plt.close("all")
    print(f"wrote {args.out}/clusters_overview.png, split_tree.png, clusters_free_energy.png, groups.csv, labels.csv")


if __name__ == "__main__":
    main()
