#!/usr/bin/env python
"""
Cluster network from the symmetry-corrected inter-cluster RMSD (Shu-Yu's suggestion): one node per
cluster, one edge per pair of clusters carrying their mean heavy-atom RMSD, drawn with networkx.

Reads a cluster_prune.py output folder (MD-tools 0.6.5.1 run):
  intercluster_rmsd.json   corrected_mean / corrected_std (A), diagonal = within-cluster spread
  merged_groups.json       cluster populations (masses)
  crystal_assignment.json  (optional) the cluster each crystal ligand is assigned to

The figure:
  nodes   area ~ population; label = cluster, population, spread (within-cluster RMSD)
  edges   colour and width = between-cluster RMSD (thick/dark = similar structures)
  layout  Kamada-Kawai with the RMSD matrix as target distances, so how far apart two nodes are
          on the page approximates their RMSD (the scale bar gives 1 A)
  many clusters: every pair is an edge, but only the minimum spanning tree (each state's chain
          of nearest neighbours) is drawn strongly; RMSD values are printed when <= 6 clusters
  gold ring = cluster the crystal is assigned to (fit.assign_reference), as in clusters_overview.png
Also writes the graph as GraphML (node: mass, within_rmsd, crystal; edge: rmsd, rmsd_std,
rmsd_by_index) for further analysis.

Needs networkx, matplotlib, numpy (env analysis-env-main).
  python cluster_network.py <clustering-dir> [<clustering-dir> ...] [--name NAME] [--out-dir DIR]
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

BEST_MATCH_COLOUR = "#e6a100"   # same gold as the frame in clusters_overview.png
NODE_COLOUR = "#86b6ef"
EDGE_CMAP = "viridis"           # dark = small RMSD = similar


def load(folder):
    with open(os.path.join(folder, "intercluster_rmsd.json")) as f:
        rmsd = json.load(f)
    with open(os.path.join(folder, "merged_groups.json")) as f:
        groups = json.load(f)
    crystal = {}
    path = os.path.join(folder, "crystal_assignment.json")
    if os.path.exists(path):
        with open(path) as f:
            for label, v in json.load(f).items():
                if v.get("best_cluster"):
                    crystal.setdefault(v["best_cluster"], []).append(label)
    return rmsd, groups, crystal


def build_graph(rmsd, groups, crystal):
    names = rmsd["clusters"]
    mean = np.asarray(rmsd["corrected_mean"], dtype=float)
    std = np.asarray(rmsd["corrected_std"], dtype=float)
    by_index = np.asarray(rmsd["by_index_mean"], dtype=float)
    g = nx.Graph()
    for i, c in enumerate(names):
        g.add_node(c, mass=float(groups["masses"][c]), within_rmsd=float(mean[i, i]),
                   crystal=",".join(crystal.get(c, [])))
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            g.add_edge(names[i], names[j], rmsd=float(mean[i, j]), rmsd_std=float(std[i, j]),
                       rmsd_by_index=float(by_index[i, j]))
    return g


def draw(g, path, title):
    nodes = list(g.nodes)
    n = len(nodes)
    rm = {(u, v): d["rmsd"] for u, v, d in g.edges(data=True)}
    dist = {u: {v: (0.0 if u == v else rm.get((u, v), rm.get((v, u)))) for v in nodes} for u in nodes}
    pos = nx.kamada_kawai_layout(g, dist=dist) if n > 2 else (
        {nodes[0]: np.array([0.0, 0.0]), nodes[1]: np.array([dist[nodes[0]][nodes[1]], 0.0])} if n == 2
        else {nodes[0]: np.array([0.0, 0.0])})
    # kamada_kawai rescales to the unit box: put it back in A, so page distance ~ RMSD
    if n >= 2:
        pts = np.array([pos[c] for c in nodes])
        page = np.array([np.linalg.norm(pts[i] - pts[j]) for i in range(n) for j in range(i + 1, n)])
        target = np.array([dist[nodes[i]][nodes[j]] for i in range(n) for j in range(i + 1, n)])
        scale = float((page * target).sum() / max((page ** 2).sum(), 1e-12))
        pos = {c: p * scale for c, p in pos.items()}

    mst = set(nx.minimum_spanning_tree(g, weight="rmsd").edges()) if n > 1 else set()
    in_mst = lambda e: e in mst or (e[1], e[0]) in mst
    values = [d["rmsd"] for _, _, d in g.edges(data=True)]
    lo, hi = (min(values), max(values)) if values else (0.0, 1.0)
    norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi if hi > lo else lo + 1e-6)
    cmap = plt.get_cmap(EDGE_CMAP)

    fig, ax = plt.subplots(figsize=(7.5 + 0.25 * n, 6.5 + 0.2 * n), layout="constrained")
    many = n > 6
    for u, v, d in g.edges(data=True):
        strong = in_mst((u, v)) or not many
        width = 1.0 + 6.0 * (hi - d["rmsd"]) / (hi - lo) if hi > lo else 3.0
        ax.plot(*zip(pos[u], pos[v]), color=cmap(norm(d["rmsd"])), lw=width if strong else 0.8,
                alpha=0.95 if strong else 0.18, zorder=1, solid_capstyle="round")
        if not many:
            mid = (pos[u] + pos[v]) / 2
            ax.text(*mid, f"{d['rmsd']:.2f}", fontsize=9, ha="center", va="center", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    masses = np.array([g.nodes[c]["mass"] for c in nodes])
    sizes = 300 + 5000 * masses / max(masses.max(), 1e-12)
    for c, s in zip(nodes, sizes):
        best = bool(g.nodes[c]["crystal"])
        ax.scatter(*pos[c], s=s, color=NODE_COLOUR, edgecolors=BEST_MATCH_COLOUR if best else "white",
                   linewidths=4 if best else 1.5, zorder=4)
        label = f"{c}\n{100 * g.nodes[c]['mass']:.1f}%\nspread {g.nodes[c]['within_rmsd']:.2f} Å"
        ax.text(*pos[c], label, ha="center", va="center", fontsize=8.5, zorder=5)
    # 1 A scale bar, so the layout can be read as a distance
    if n >= 2:
        xs = np.array([p[0] for p in pos.values()]); ys = np.array([p[1] for p in pos.values()])
        x0, y0 = xs.min(), ys.min() - 0.12 * max(np.ptp(ys), np.ptp(xs), 1.0)
        ax.plot([x0, x0 + 1.0], [y0, y0], color="#333", lw=2)
        ax.text(x0 + 0.5, y0, "1 Å", ha="center", va="bottom", fontsize=9, color="#333")
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    if values:
        fig.colorbar(sm, ax=ax, shrink=0.7, label="between-cluster RMSD, symmetry-corrected (Å)")
    ax.set_aspect("equal")
    ax.margins(0.18)
    ax.axis("off")
    crystals = [f"{lab} → {c}" for c in nodes for lab in g.nodes[c]["crystal"].split(",") if lab]
    sub = ("node area = population, spread = within-cluster RMSD; page distance ≈ RMSD"
           + ("; strong edges = minimum spanning tree" if many else "")
           + ("\ngold ring = crystal assigned: " + ", ".join(crystals) if crystals else ""))
    ax.set_title(f"{title}\n{sub}", fontsize=10)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("folders", nargs="+", help="cluster_prune.py output folder(s)")
    ap.add_argument("--name", default=None, help="title (default: the folder name)")
    ap.add_argument("--out-dir", default=None, help="write here instead of into each folder")
    args = ap.parse_args()
    for folder in args.folders:
        folder = os.path.abspath(folder)
        if not os.path.exists(os.path.join(folder, "intercluster_rmsd.json")):
            print(f"skip {folder}: no intercluster_rmsd.json (no cluster, or an older run)")
            continue
        rmsd, groups, crystal = load(folder)
        g = build_graph(rmsd, groups, crystal)
        name = args.name or os.path.basename(folder.rstrip("/"))
        out = args.out_dir or folder
        os.makedirs(out, exist_ok=True)
        stem = os.path.join(out, f"{name}_cluster_network" if args.out_dir else "cluster_network")
        nx.write_graphml(g, stem + ".graphml")
        draw(g, stem + ".png", f"{name}: cluster network\n{g.number_of_nodes()} clusters, "
                               f"heavy-atom RMSD over {rmsd['heavy_atoms']} atoms, {rmsd['n_operations']} symmetry operations")
        print(f"wrote {stem}.png/.graphml  ({g.number_of_nodes()} nodes, {g.number_of_edges()} edges)")


if __name__ == "__main__":
    sys.exit(main())
