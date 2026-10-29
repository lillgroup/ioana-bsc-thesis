#!/usr/bin/env python
"""
Recursive symmetry-first clustering: cluster, and re-cluster every cluster whose internal
(within-cluster) symmetry-corrected heavy-atom RMSD is above --max-rmsd, on that cluster's frames
alone, until no final cluster is above the limit or a cluster cannot be split further.

Each level is the same fit as cluster_prune.py (MD-tools 0.6.5.1 t_hdbscan via build_fit, default
allow_single_cluster=False), with a stricter vote if asked (--k, --min-vote-fraction: 19/20 is
k=20, 0.95; 39/40 is k=40, 0.975). The internal RMSD is intercluster_rmsd's diagonal (mean RMSD
between frames of the same cluster, minimised over the symmetry operations).

A cluster is final ("leaf") when its internal RMSD <= --max-rmsd, or when re-clustering it gives
no cluster (one state) or only one cluster, or --max-depth is reached. Frames that a re-clustering
leaves unassigned (noise / ambiguous) are unassigned in the final result. Names follow the tree:
A, then A.A, A.B, ... (A.B = sub-cluster B of cluster A).

The mass floor is a fraction of the frames being clustered; for a subset it is raised to at least
--k / n frames, because a cluster needs about k training frames to win the vote (a smaller one is
emptied by it and MD-tools' alignment fails).

Outputs in --out: tree.csv (every node: depth, frames, mass of the whole ensemble, internal RMSD,
leaf or split), labels.csv (final leaf per frame), recursive_overview.png (every leaf in one
figure, each in the labelling of its own fit), recursive.log via stdout.
Needs env analysis-env-main.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import mdtraj as md
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_prune import load_torsions, build_fit, draw_overview, crystal_rows_to_refs  # noqa: E402
from md_tools.analysis import (group_name, enumerate_symmetry, heavy_atoms,  # noqa: E402
                               operations_on_atoms, intercluster_rmsd, reference_assignments)


def fit_subset(tors, cmd, mol, torsions, idx, args):
    """Cluster the frames idx; returns (fit, aligned angles of the selected torsions) or None."""
    n = len(idx)
    floor = max(args.mass_floor, (args.k or 20) / n)
    for _ in range(4):                       # if the vote empties a cluster, raise the floor
        try:
            fit, mapping, _ = build_fit(tors[idx], cmd[idx], mol, torsions, mass_floor=floor,
                                        max_matches=args.max_matches, k=args.k,
                                        min_vote_fraction=args.min_vote_fraction)
            al = fit.align_all(per_cluster=True)
            return fit, al["selected"], mapping, floor
        except RuntimeError as e:
            if "no draw produced a cluster" in str(e):
                return None
            raise
        except ValueError as e:
            if "has no members" not in str(e):
                raise
            floor *= 2
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--sdf", required=True)
    ap.add_argument("--torsions-file", action="append", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--mass-floor", type=float, default=0.005)
    ap.add_argument("--k", type=int, default=None)
    ap.add_argument("--min-vote-fraction", type=float, default=None)
    ap.add_argument("--max-matches", type=int, default=None)
    ap.add_argument("--max-rmsd", type=float, default=2.0, help="split clusters above this internal RMSD (A)")
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--rmsd-max-per-cluster", type=int, default=1500)
    ap.add_argument("--ref-sdf", action="append", default=None)
    ap.add_argument("--name", default="ligand")
    ap.add_argument("--log-floor", type=float, default=None,
                    help="also draw recursive_overview_log.png: fraction of frames per bin on a log "
                         "axis cut at this value (e.g. 1e-5)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    torsions = load_torsions(args.torsions_file)
    names, quads = list(torsions), list(torsions.values())
    cmd = md.load(args.traj, top=args.top, stride=args.stride)
    tors = md.compute_dihedrals(cmd, quads)
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    n_all = cmd.n_frames
    heavy = heavy_atoms(mol)
    sym_kw = {} if args.max_matches is None else {"max_matches": int(args.max_matches)}
    operations = operations_on_atoms(enumerate_symmetry(mol, **sym_kw), heavy)
    print(f"{n_all} frames, {len(names)} torsions; vote k={args.k or 20}, "
          f"min_vote_fraction={args.min_vote_fraction or 0.90}; split while internal RMSD > {args.max_rmsd} A")

    final = np.full(n_all, "unassigned", dtype=object)
    leaves, nodes = {}, []
    descriptor = None
    stack = [("", np.arange(n_all), 0)]          # (node name, frame indices, depth)
    while stack:
        node, idx, depth = stack.pop(0)
        res = fit_subset(tors, cmd, mol, torsions, idx, args)
        label = node or "(all)"
        if res is None:
            print(f"{'  ' * depth}{label}: {len(idx)} frames -> no cluster (one state); kept as a leaf")
            if node:
                leaves[node]["note"] = "cannot split (one state)"
            continue
        fit, sel, mapping, floor = res
        descriptor = descriptor or fit.descriptor_
        X = cmd.xyz[idx][:, np.asarray(mapping)[heavy], :] * 10.0
        r = intercluster_rmsd(X, fit.labels_, operations, max_per_cluster=args.rmsd_max_per_cluster,
                              workers=os.cpu_count() or 1)
        within = {c: float(r["mean"][i, i]) for i, c in enumerate(r["clusters"])}
        print(f"{'  ' * depth}{label}: {len(idx)} frames, mass floor {floor:.4f} -> {fit.n_clusters_} clusters, "
              f"unassigned {100 * fit.unassigned_mass_:.1f}%")
        if node and fit.n_clusters_ < 2:
            print(f"{'  ' * depth}  only one sub-cluster: {label} kept as a leaf")
            leaves[node]["note"] = "cannot split (one sub-cluster)"
            continue
        if node:                                     # node is split: its frames move to the children
            final[idx] = "unassigned"
            leaves.pop(node)
        for c in range(fit.n_clusters_):
            child = (node + "." if node else "") + group_name(c)
            members = idx[fit.labels_ == c]
            rmsd_c = within[c]
            nodes.append({"cluster": child, "parent": node or "", "depth": depth + 1,
                          "frames": len(members), "mass_of_ensemble": len(members) / n_all,
                          "internal_rmsd_A": round(rmsd_c, 3)})
            final[members] = child
            leaves[child] = {"frames": members, "rmsd": rmsd_c, "angles": sel[fit.labels_ == c], "note": ""}
            split = rmsd_c > args.max_rmsd and depth + 1 < args.max_depth
            print(f"{'  ' * depth}  {child}: {100 * len(members) / n_all:5.1f}% of ensemble, internal RMSD "
                  f"{rmsd_c:.2f} A" + ("  -> split" if split else ""))
            if split:
                stack.append((child, members, depth + 1))

    tree = pd.DataFrame(nodes)
    tree["leaf"] = tree["cluster"].isin(list(leaves))
    tree["note"] = tree["cluster"].map(lambda c: leaves.get(c, {}).get("note", "") if c in leaves else "split")
    tree.to_csv(os.path.join(args.out, "tree.csv"), index=False)
    pd.DataFrame({"frame": np.arange(n_all), "time_ps": cmd.time, "cluster": final}).to_csv(
        os.path.join(args.out, "labels.csv"), index=False)

    order = sorted(leaves, key=lambda c: -len(leaves[c]["frames"]))
    print("\nfinal clusters:")
    for c in order:
        v = leaves[c]
        print(f"  {c:10s} {100 * len(v['frames']) / n_all:5.1f}%  internal RMSD {v['rmsd']:.2f} A  {v['note']}")
    unassigned = float(np.mean(final == "unassigned"))
    above = [c for c in order if leaves[c]["rmsd"] > args.max_rmsd]
    print(f"unassigned {100 * unassigned:.1f}%; leaves above {args.max_rmsd} A: {above or 'none'}")

    refs = {}
    for path in (args.ref_sdf or []):
        rows = reference_assignments(path, mol, descriptor.quads)
        refs[os.path.splitext(os.path.basename(path))[0]] = crystal_rows_to_refs(rows, len(names))
    aligned = {c: leaves[c]["angles"] for c in order}
    masses = [len(leaves[c]["frames"]) / n_all for c in order]
    labels = [f"{c}\n({100 * m:.1f}%)\nRMSD {leaves[c]['rmsd']:.2f} Å" + ("\n" + leaves[c]["note"] if leaves[c]["note"] else "")
              for c, m in zip(order, masses)]
    vote = f"{int(round((args.min_vote_fraction or 0.90) * (args.k or 20)))}/{args.k or 20} vote"
    title = (f"{args.name} — recursive clustering ({vote}, split while internal RMSD > "
             f"{args.max_rmsd} Å): {len(order)} final clusters, unassigned {100 * unassigned:.1f}%\n"
             "each cluster row in the labelling of its own fit; black lines = crystal, every assignment")
    for fname, floor in (("recursive_overview.png", None), ("recursive_overview_log.png", args.log_floor)):
        if fname.endswith("_log.png") and not floor:
            continue
        draw_overview(os.path.join(args.out, fname), mol=mol, names=names, quads=quads,
                      tors=tors, aligned=aligned, masses=masses, noise=unassigned, refs=refs, row_labels=labels,
                      top_label="all frames,\nas stored", log_floor=floor,
                      title=title + (f"; y = fraction of the row's frames per bin, log, cut at {floor:g}" if floor else ""))
    print(f"wrote {args.out}/ (tree.csv, labels.csv, recursive_overview[_log].png)")


if __name__ == "__main__":
    main()
