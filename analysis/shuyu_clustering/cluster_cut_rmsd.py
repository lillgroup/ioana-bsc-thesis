#!/usr/bin/env python
"""
Inter-cluster RMSD and cluster network for a TORSION PARTITIONING (barrier-cut) result: the same two
figures cluster_prune.py + cluster_network.py make for the symmetry-first HDBSCAN clusters, on the
clusters of a finished cluster_cut.py run (its labels.csv and summary.json; nothing is re-clustered).

Barrier cuts label every frame and can give hundreds of small cells, so only clusters holding
>= --min-mass of all frames (default 1 %, as the clusters drawn in clusters_overview.png) become rows
of the matrix / nodes of the network; the frames of smaller clusters are left out, and the title says
how much that is. RMSD: symmetry-corrected heavy-atom RMSD (md_tools intercluster_rmsd, at most
--rmsd-max-per-cluster frames per cluster, as cluster_prune.py), diagonal = within-cluster spread.

Writes into --cut-dir: intercluster_rmsd.csv/.json/.png, merged_groups.json (populations) and
crystal_assignment.json (cluster_cut's crystal cluster, for the gold ring), cluster_network.png/.graphml.

  python cluster_cut_rmsd.py --cut-dir <cluster_cut --out> --traj ... --top ... --sdf ... [--max-matches N] [--name NAME]
"""
import argparse
import json
import os
import sys

import mdtraj as md
import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_prune import write_intercluster_rmsd  # noqa: E402
from cluster_network import load, build_graph, draw  # noqa: E402
from md_tools.analysis import (verify_atom_mapping, enumerate_symmetry, heavy_atoms,  # noqa: E402
                               operations_on_atoms)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cut-dir", required=True, help="cluster_cut.py output folder")
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True, help="topology WITH bonds")
    ap.add_argument("--sdf", required=True)
    ap.add_argument("--max-matches", type=int, default=None)
    ap.add_argument("--min-mass", type=float, default=0.01, help="smallest cluster shown, fraction of all frames")
    ap.add_argument("--rmsd-max-per-cluster", type=int, default=1500, help="0 = every frame")
    ap.add_argument("--rmsd-workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--name", default=None)
    args = ap.parse_args()
    out = os.path.abspath(args.cut_dir)
    P = lambda f: os.path.join(out, f)
    name = args.name or os.path.basename(out)

    summary = json.load(open(P("summary.json")))
    masses = summary["cluster_mass"]                                 # name -> fraction, by mass
    shown = [c for c, m in masses.items() if m >= args.min_mass]
    index = {c: i for i, c in enumerate(shown)}
    labels_csv = pd.read_csv(P("labels.csv"), keep_default_na=False, dtype={"cluster": str})
    labels = np.array([index.get(c, -1) for c in labels_csv["cluster"]], dtype=np.int64)

    cmd = md.load(args.traj, top=args.top, stride=int(summary["stride"]))
    if cmd.n_frames != len(labels):
        sys.exit(f"{cmd.n_frames} frames, labels.csv has {len(labels)}")
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    mapping, _ = verify_atom_mapping(mol, [(a.index, b.index) for a, b in cmd.topology.bonds],
                                     [a.element.symbol for a in cmd.topology.atoms])
    heavy = heavy_atoms(mol)
    sym_kw = {} if args.max_matches is None else {"max_matches": int(args.max_matches)}
    operations = operations_on_atoms(enumerate_symmetry(mol, **sym_kw), heavy)
    X = cmd.xyz[:, np.asarray(mapping)[heavy], :] * 10.0             # nm -> A, heavy atoms in the molecule's order

    rest = 1.0 - sum(masses[c] for c in shown)
    n_small = len(masses) - len(shown)
    left_out = f"{n_small} clusters < {100 * args.min_mass:g} % left out ({100 * rest:.1f} % of frames)" if n_small else ""
    print(f"{name}: {len(shown)} of {len(masses)} barrier-cut clusters >= {100 * args.min_mass:g} %; {left_out}")
    write_intercluster_rmsd(P, f"{name} (barrier cuts{', ' + left_out if left_out else ''}): ", X, labels,
                            lambda c: shown[c], operations, args.rmsd_max_per_cluster or None, args.rmsd_workers)

    json.dump({"masses": {c: masses[c] for c in shown}}, open(P("merged_groups.json"), "w"), indent=2)
    json.dump({lab: {"best_cluster": v["cluster"] if v["cluster"] in index else None, "note": v.get("note", "")}
               for lab, v in summary.get("crystal_assignment", {}).items()},
              open(P("crystal_assignment.json"), "w"), indent=2)
    rmsd, groups, crystal = load(out)
    g = build_graph(rmsd, groups, crystal)
    nx.write_graphml(g, P("cluster_network.graphml"))
    draw(g, P("cluster_network.png"),
         f"{name}: cluster network, barrier cuts\n{g.number_of_nodes()} clusters >= {100 * args.min_mass:g} %"
         f"{' (' + left_out + ')' if left_out else ''}, heavy-atom RMSD over {rmsd['heavy_atoms']} atoms, "
         f"{rmsd['n_operations']} symmetry operations")
    print(f"wrote {out}/intercluster_rmsd.png, cluster_network.png")


if __name__ == "__main__":
    main()
