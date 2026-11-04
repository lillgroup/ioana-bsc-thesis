"""The sampled torsion distributions of one system, with the crystal value: the top row of the
clustering figures on its own (no clusters).

One panel per torsion t0.. (prepare.py's numbering, the same as every other figure of the compound),
REST2 state 0, every frame, written in ONE labelling: each frame relabelled by the molecule's symmetry
to the labelling nearest frame 0 (descriptor.align_to_reference, as cut_clustering's aligned_angles),
so a symmetric group (a carboxylate, a ring flip) shows one peak, not two. Crystal lines: the crystal
ligand in the labelling nearest that same frame.

  python -I sampled_torsions.py --name folic-acid --system 3QL3 [--out PREFIX]
default output: results/<name>/<version>/<name>_<system>_sampled_torsions.png
"""
import argparse
import io
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "shuyu_clustering"))
from cluster_prune import (load_torsions, molecule_png, torsion_colours, draw_reference_lines,  # noqa: E402
                           read_rest2_config, simulated_ns)
from md_tools.analysis import t_hdbscan, TorsionDefinitions, verify_atom_mapping, reference_assignments  # noqa: E402

WORK = "/scicore/home/lill0000/popa0002/thesis/analysis-report"
REPO = os.path.dirname(os.path.dirname(HERE))
KINDS = ("torsions", "torsions_rigid", "torsions_polarH")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True, help="compound (compounds.csv)")
    ap.add_argument("--system", required=True, help="water, or a complex label (e.g. 3QL3)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    man = json.load(open(f"{WORK}/{args.name}/manifest.json"))
    sys_ = man["water"] if args.system == "water" else next(c for c in man["complexes"] if c["label"] == args.system)
    out = args.out or os.path.join(REPO, "results", args.name, man["version"],
                                   f"{args.name}_{args.system}_sampled_torsions")

    tors = load_torsions([f"{sys_['torsions']}/{k}.txt" for k in KINDS if os.path.getsize(f"{sys_['torsions']}/{k}.txt")])
    names, quads = list(tors), list(tors.values())
    t = md.load(sys_["traj"], top=sys_["top"])
    mol = Chem.SDMolSupplier(sys_["sdf"], removeHs=False)[0]
    mapping, _ = verify_atom_mapping(mol, [(a.index, b.index) for a, b in t.topology.bonds],
                                     [a.element.symbol for a in t.topology.atoms])
    defs = TorsionDefinitions(quads=tuple(quads), names=tuple(names), indexing="zero-based", units="radian")
    extra = {"max_matches": int(man["max_matches"])} if man["max_matches"] else {}
    fit = t_hdbscan(md.compute_dihedrals(t, quads), molecule=mol, torsion_definitions=defs, coordinates=t,
                    atom_mapping=mapping, **extra)                       # not fitted: symmetry only
    d = fit.descriptor_
    C = np.asarray(fit.closed_theta_, float)
    K = len(names)
    wrap = lambda x: np.angle(np.exp(1j * x))
    angles = wrap(d.align_to_reference(C, C[0])["angles"][:, :K])
    refs = {}
    for path in sys_["refs"]:
        near = d.nearest_labelling(reference_assignments(path, mol, d.quads), C[0])["angles"]
        refs[os.path.splitext(os.path.basename(path))[0]] = [[float(wrap(a))] if np.isfinite(a) else [] for a in near[:K]]

    colours = torsion_colours(K)
    ncol = K if K <= 6 else -(-K // 2)              # one row, or two for many torsions (readable panels)
    nrow = -(-K // ncol)
    panel_w, panel_h, mol_w = 2.5, 2.4, 5.5
    fig = plt.figure(figsize=(mol_w + ncol * panel_w, max(nrow * panel_h, 3.6) + 0.6), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=[mol_w, ncol * panel_w])
    mol_ax = fig.add_subplot(outer[0])
    mol_ax.imshow(mpimg.imread(io.BytesIO(molecule_png(mol, quads, colours, size=(1100, 1100))), format="png"))
    mol_ax.axis("off")
    mol_ax.set_title("t<k> marks the central bond of torsion k", fontsize=10)
    grid = outer[1].subgridspec(nrow, ncol)
    edges = np.linspace(-np.pi, np.pi, 73)
    for k in range(K):
        ax = fig.add_subplot(grid[divmod(k, ncol)])
        ax.hist(angles[:, k], bins=edges, density=True, color=colours[k], alpha=0.85)
        draw_reference_lines(ax, refs, k, with_label=(k == 0))
        ax.set_xlim(-np.pi, np.pi)
        ax.set_xticks([-np.pi, 0, np.pi], ["-π", "0", "π"])
        ax.set_yticks([])
        kind = names[k].rsplit("_", 1)[0].replace("torsions_", "").replace("torsions", "rotatable")
        ax.set_title(f"t{k}\n{kind}", color=colours[k], fontsize=9)
        if k == 0 and refs:
            ax.legend(fontsize=7, loc="upper left", frameon=False)
    n_rep, tau = read_rest2_config(sys_.get("rest2_config"))
    label = "water" if args.system == "water" else f"{args.system} complex"
    fig.suptitle(f"{args.name} {label}: sampled torsion distributions\nREST2 state 0, {simulated_ns(t):g} ns, "
                 f"{t.n_frames} frames, one labelling; black line = crystal (bioactive) conformation", fontsize=11)
    fig.savefig(out + ".png", dpi=170)
    print(f"{args.name} {args.system} ({os.path.relpath(sys_['run'], '/scicore/home/lill0000/popa0002/thesis')}): "
          f"{t.n_frames} frames, {K} torsions -> {out}.png")


if __name__ == "__main__":
    main()
