#!/usr/bin/env python
"""
Shu-Yu's SYMMETRY-FIRST torsion clustering (MD-tools 0.6.5.1, main 2e9e205; tutorials
docs/tutorial/paracetamol/analysis/symmetry-alignment.ipynb and rgdfv/analysis/clustering.ipynb),
run on YOUR trajectory, with this thesis' figures (clusters_overview, 2-D torsion
plots, representatives, crystal assignment, inter-cluster RMSD).

Four parts, all MD-tools calls as in the tutorials:
  1. Clustering   t_hdbscan(tors, molecule, torsion_definitions, coordinates, atom_mapping):
                  the molecular symmetry is enumerated FIRST and built into the distance
                  (min over relabellings), then HDBSCAN on that dense matrix, then every frame
                  classified by an 18-of-20 neighbour vote. No merge step: a cluster is final.
                  Labels: -2 density noise, -1 ambiguous vote.
                  Since 0.6.5.1 the molecule is CANONICALISED first (torsional folding): resonance-
                  equivalent atoms (the two O of a carboxylate, the N of a guanidinium) are one
                  atom, so a torsion to one O is folded onto the torsion to the other. Default in
                  t_hdbscan, nothing to pass. allow_single_cluster now defaults to False (sklearn's).
     Gauge        --gauge <torsion> (MD-tools 0.6.5.2, tutorial section 10): for a symmetry that moves
                  several torsions together (a ring flip shifts both ring torsions by pi), the figures
                  write that ONE torsion in a half-turn window centred on the density peak and the others
                  follow it: align_all(per_cluster=True, gauge=(torsion, pi)). Drawing only; the
                  clustering distance is the minimum over every operation and is not touched.
  4. Crystal      reference_assignments(crystal, molecule, descriptor.quads) + fit.assign_reference():
     + RMSD       per cluster the crystal labelling nearest its density peak (one line per cluster,
                  not one per assignment) and which cluster the crystal resembles, or none when the
                  crystal's columns cannot tell two clusters apart. intercluster_rmsd(): mean +- std
                  heavy-atom RMSD between and within clusters, minimised over the symmetry
                  operations (and, for comparison, by atom index).
  2. Alignment    fit.align_all(per_cluster=True): every frame written in one labelling for the
                  figures (each cluster aligned to its own density peak, each peak to the density
                  peak of cluster A). Changes no label and no population.
  3. Summary      fit.representatives() (highest vote, then symmetry-aware medoid), circular
                  statistics on the aligned torsions, figures.

The previous cluster-then-merge version (t_hdbscan + TorsionSymmetry) is in git history.
Needs MD-tools main with the symmetry-first analysis (env analysis-env-main), scikit-learn, rdkit, pandas.

Example (paracetamol, state-0 solute trajectory):
  python para_clustering.py \
      --traj ~/PARA-REST2-long/REST2-run1/solute_state0_prod1.nc \
      --top  ~/PARA-REST2-long/build/built.solute.pdb \
      --sdf  ~/PARA-REST2-long/build/PAR.sdf \
      --out  ~/PARA-REST2-long/tutorial-clustering

Torsions (pick one way; indices must match the trajectory's atom numbering):
  --torsions-file torsions.txt         from select_torsions.py (any molecule, e.g. folic acid);
                                       add torsions_rigid.txt as well to include e.g. the amide
  --torsion name:i,j,k,l  (repeatable)  explicit, e.g. polar-H torsions such as C-O-H of a COOH
  neither                               the tutorial's paracetamol set (omega, aryl, phenol)
  --use-torsions 1,2,3,4                keep only these (0-based, in the order listed)

Folic acid example (solute trajectory, ligand indices):
  python para_clustering.py --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb \
      --sdf build/FOL.sdf --torsions-file REST2-run1/torsions.txt --stride 5 --out fol-clustering

Labelling and the bioactive reference (all optional):
  --name 4ibi                     compound name for the figure titles
  --rest2-config REST2.config     number of replicas and tau_max for the overview title
  --ref-sdf crystal/4ibi_1D8_A.sdf  crystal ligand (repeatable, e.g. one per chain); its torsions
                                  are drawn as vertical lines. Matched to --sdf on heavy atoms by
                                  element + connectivity only (reference_assignments), so torsions
                                  through a hydrogen get no line. Top row (all frames, as stored):
                                  every assignment; cluster rows: the assignment nearest that
                                  cluster's density peak (fit.assign_reference), as in the tutorial.
  -> crystal_assignment.csv/.json which cluster each crystal resembles, per cluster the distance to
                                  its density peak and nearest frame, and the peak percentile
  -> intercluster_rmsd.csv/.json/.png  symmetry-corrected heavy-atom RMSD between/within clusters
  -> clusters_overview_unaligned.png  the same layout with each cluster's torsions AS STORED (not
                                  relabelled): one conformer in two labellings shows both images
  -> clusters_overview.png        molecule on the left; to its right one column per torsion: top row
                                  the full distribution, then one row per cluster; many
                                  torsions wrap into 2-3 bands, many groups into 2-3 side-by-side
                                  blocks, to keep the figure ~16:9
"""
import argparse
import importlib.metadata
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md

try:
    import md_tools
    from md_tools.analysis import (t_hdbscan, TorsionDefinitions, verify_atom_mapping, group_name,
                                   strip_nonpolar_hydrogens, draw_representative_structures,
                                   reference_assignments, enumerate_symmetry, heavy_atoms,
                                   operations_on_atoms, intercluster_rmsd)
    from md_tools.analysis import THDBSCAN
    if not hasattr(THDBSCAN, "assign_reference"):
        raise ImportError("this md_tools has no THDBSCAN.assign_reference (needs 0.6.5.1)")
except ImportError as e:
    sys.exit(f"md_tools.analysis not available ({e}).\n"
             "Needs Shu-Yu's MD-tools 0.6.5.1 (main 2e9e205): env analysis-env-main.")
from rdkit import Chem
import pandas as pd

DEFAULT_TORSIONS = ["omega:0,1,3,4", "aryl:1,3,4,5", "phenol:6,7,8,17"]

# One colour per TORSION, all distinct (tab10 without its grey, then tab20's light shades).
# MD-tools' own palette has 5 entries, so with 7 torsions two pairs shared a colour.
TORSION_PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                   "#e377c2", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a",
                   "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5"]
REF_STYLES = ["-", "--", ":", "-."]


def circular_mean(values):
    return float(np.arctan2(np.sin(values).mean(), np.cos(values).mean()))


def circular_std(values):
    resultant = np.hypot(np.sin(values).mean(), np.cos(values).mean())
    return float(np.sqrt(-2.0 * np.log(max(resultant, 1e-12))))


def torsion_colours(n):
    if n > len(TORSION_PALETTE):
        sys.exit(f"{n} torsions but only {len(TORSION_PALETTE)} distinct colours defined")
    return TORSION_PALETTE[:n]


def load_torsions(files=None, specs=None, use=None):
    """name -> (i, j, k, l), from select_torsions.py files and/or name:i,j,k,l specs."""
    torsions = {}
    for path in (files or []):
        for line in open(path):
            if line.strip().startswith("#") or not line.strip():
                continue
            q = tuple(int(x) for x in line.split("#")[0].split()[:4])
            label = os.path.splitext(os.path.basename(path))[0]
            torsions[f"{label}_{len(torsions)}"] = q
    for spec in (specs or ([] if torsions else DEFAULT_TORSIONS)):
        name, idx = spec.split(":")
        torsions[name] = tuple(int(x) for x in idx.split(","))
    if use:
        keep_pos = [int(x) for x in use.split(",")]
        items = list(torsions.items())
        torsions = dict(items[i] for i in keep_pos)
    if not torsions:
        sys.exit("no torsions given")
    return torsions


def read_rest2_config(path):
    """(number_of_replicas, tau_max) from a REST2.config, or (None, None)."""
    if not path:
        return None, None
    import yaml
    with open(path) as f:
        rest2 = (yaml.safe_load(f) or {}).get("rest2", {})
    return rest2.get("number_of_replicas"), rest2.get("tau_max")


def simulated_ns(traj, stride=1):
    """Simulated time covered by a trajectory, from its own time stamps (ps).

    A strided load drops up to stride-1 frames at the end (60 ns at stride 5 ends at
    59.992 ns), so they are added back from the unstrided frame interval."""
    t = np.asarray(traj.time, dtype=float)
    if t.size < 2:
        return None
    dt = (t[1] - t[0]) / stride
    return float(t[-1] + (stride - 1) * dt) / 1000.0


def _heavy_graph(mol):
    """Heavy-atom copy with every bond single and no aromaticity: matches on element +
    connectivity only, so it does not depend on how the crystal file assigns bond orders."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rw = Chem.RWMol()
    new = {}
    for i in heavy:
        atom = Chem.Atom(mol.GetAtomWithIdx(i).GetAtomicNum())
        atom.SetNoImplicit(True)
        new[i] = rw.AddAtom(atom)
    for b in mol.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in new and c in new:
            rw.AddBond(new[a], new[c], Chem.BondType.SINGLE)
    g = rw.GetMol()
    g.UpdatePropertyCache(strict=False)
    return g, heavy


def reference_torsions(ref_sdf, sim_mol, quads):
    """Torsions of a crystal ligand, in the simulation's atom numbering.

    Returns one list per torsion: the distinct values (rad) over every graph match, or [] when
    the torsion goes through a hydrogen (crystal structures carry none)."""
    ref = Chem.MolFromMolFile(ref_sdf, removeHs=False, sanitize=False)
    if ref is None:
        sys.exit(f"RDKit could not read {ref_sdf}")
    sim_g, sim_heavy = _heavy_graph(sim_mol)
    ref_g, ref_heavy = _heavy_graph(ref)
    if sim_g.GetNumAtoms() != ref_g.GetNumAtoms():
        sys.exit(f"{ref_sdf}: {ref_g.GetNumAtoms()} heavy atoms, simulation ligand has "
                 f"{sim_g.GetNumAtoms()}")
    matches = sim_g.GetSubstructMatches(ref_g, uniquify=False, useChirality=False,
                                        maxMatches=100000)
    if not matches:
        sys.exit(f"{ref_sdf}: heavy-atom graph does not match the simulation ligand")
    xyz = ref.GetConformer().GetPositions()
    sim_pos = {atom: k for k, atom in enumerate(sim_heavy)}

    # The graph match ignores chirality, so for e.g. two identical aryls on two stereocentres it
    # also returns the mirror-image assignment. Keep only matches under which every stereocentre
    # of the simulated ligand has the same handedness in the crystal.
    sim_xyz = sim_mol.GetConformer().GetPositions()
    probe = Chem.Mol(sim_mol)
    Chem.AssignStereochemistryFrom3D(probe)
    centres = []
    for c, _ in Chem.FindMolChiralCenters(probe, includeUnassigned=False):
        nbrs = [n.GetIdx() for n in probe.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() > 1]
        if len(nbrs) >= 3:
            centres.append((c, nbrs[:3], _handedness(sim_xyz, c, nbrs[:3])))
    if centres:
        def keeps(m):
            inv = {s: r for r, s in enumerate(m)}
            pick = lambda a: ref_heavy[inv[sim_pos[a]]]
            return all(_handedness(xyz, pick(c), [pick(n) for n in nb]) == sign
                       for c, nb, sign in centres)
        kept = [m for m in matches if keeps(m)]
        if not kept:
            print(f"  WARNING {ref_sdf}: crystal ligand is not the simulated stereoisomer; "
                  f"reference torsions ignore chirality")
        else:
            matches = kept
    out = []
    for q in quads:
        if any(a not in sim_pos for a in q):
            out.append([])
            continue
        values = []
        for m in matches:            # m[ref_graph_idx] = sim_graph_idx
            inv = {s: r for r, s in enumerate(m)}
            p = xyz[[ref_heavy[inv[sim_pos[a]]] for a in q]]
            v = _dihedral(p)
            if all(abs(np.angle(np.exp(1j * (v - w)))) > np.deg2rad(2) for w in values):
                values.append(v)
        out.append(values)
    return out


def _handedness(xyz, centre, nbrs):
    """Sign of the volume spanned by three neighbours around a stereocentre."""
    v = xyz[nbrs] - xyz[centre]
    return bool(np.linalg.det(v) > 0)


def _dihedral(p):
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def draw_reference_lines(axis, refs, column, with_label=False, colour="black"):
    """Vertical lines at the crystal value(s) of torsion `column`, one line style per file. A value
    at the seam (|v| > 3.0 rad) shows at both ends, as in Shu-Yu's tutorial."""
    for k, (label, values) in enumerate(refs.items()):
        first = with_label
        for v in values[column]:
            for x in ((v, v - 2 * np.pi * np.sign(v)) if abs(v) > 3.0 else (v,)):
                axis.axvline(x, color=colour, lw=1.3, ls=REF_STYLES[k % len(REF_STYLES)],
                             label=(f"crystal {label}" if first else None))
                first = False


def crystal_rows_to_refs(rows, n_selected):
    """reference_assignments rows -> {torsion column: distinct finite values} for the selected
    torsions (every assignment of the crystal, as the tutorial's 'as stored' row draws them)."""
    out = []
    for t in range(n_selected):
        values = []
        for v in rows[:, t]:
            if np.isfinite(v) and all(abs(np.angle(np.exp(1j * (v - w)))) > np.deg2rad(2) for w in values):
                values.append(float(v))
        out.append(values)
    return out


def molecule_png(mol, quads, colours, size=(900, 900)):
    """The ligand (non-polar H stripped) with each torsion's three bonds highlighted in its own
    colour. Bonds or atoms shared by several torsions are drawn in all their colours, so no
    torsion hides another; the central bond carries the torsion number."""
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    keep_h = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1
              and any(n.GetAtomicNum() not in (1, 6) for n in a.GetNeighbors())}
    keep_h |= {a for q in quads for a in q if mol.GetAtomWithIdx(a).GetAtomicNum() == 1}
    rw = Chem.RWMol(mol)
    for a in sorted((a.GetIdx() for a in mol.GetAtoms()
                     if a.GetAtomicNum() == 1 and a.GetIdx() not in keep_h), reverse=True):
        rw.RemoveAtom(a)
    draw = rw.GetMol()
    old_to_new, n = {}, 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 1 or a.GetIdx() in keep_h:
            old_to_new[a.GetIdx()] = n
            n += 1
    draw.UpdatePropertyCache(strict=False)
    AllChem.Compute2DCoords(draw)

    from matplotlib.colors import to_rgb
    atom_hl, bond_hl = {}, {}
    for k, (q, colour) in enumerate(zip(quads, colours)):
        rgb = to_rgb(colour)
        q = [old_to_new[a] for a in q]
        for a in q:
            atom_hl.setdefault(a, []).append(rgb)
        for x, y in zip(q, q[1:]):
            bond = draw.GetBondBetweenAtoms(x, y)
            if bond is not None:
                bond_hl.setdefault(bond.GetIdx(), []).append(rgb)
        centre = draw.GetBondBetweenAtoms(q[1], q[2])
        if centre is not None:
            note = centre.GetProp("bondNote") + "," if centre.HasProp("bondNote") else ""
            centre.SetProp("bondNote", f"{note}t{k}")
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.annotationFontScale = 0.9
    opts.fillHighlights = False
    opts.continuousHighlight = False
    drawer.DrawMoleculeWithHighlights(draw, "", atom_hl, bond_hl, {}, {})
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


BEST_MATCH_COLOUR = "#e6a100"   # frame around the group the crystal is assigned to


def draw_overview(path, *, mol, names, quads, tors, aligned, masses, noise, refs, title,
                  row_labels=None, cluster_refs=None, top_refs=None, top_label="all frames",
                  best=None, log_floor=None):
    """One figure for the whole clustering: molecule on the left; to its right one column per
    torsion and one row per population -- the full distribution first, then each merged group --
    all in the torsion's colour, crystal lines on top. The grid is wrapped to keep the whole
    figure near 16:9: many torsions with few groups -> the torsions wrap into 1-3 bands stacked
    vertically (each band repeats the population rows); many groups with few torsions -> the
    populations wrap into 1-3 blocks side by side (each block repeats the torsion columns).
    `cluster_refs` (group -> refs) gives each cluster row its own crystal lines (the labelling
    nearest that cluster's density peak); without it every row shows `refs`. `tors` is what the
    top row draws (stored or aligned angles) and `top_refs` its crystal lines (default `refs`).
    `best` (group -> crystal labels) frames the group each crystal is assigned to
    (fit.assign_reference's best_cluster) and names the crystal in that row's label."""
    import io
    import matplotlib.image as mpimg

    colours = torsion_colours(len(names))
    groups = list(aligned)
    nrow = 1 + len(groups)
    ntor = len(names)
    gap, target = 0.5, 16 / 9

    def size(nband, nblock, panel_h, panel_w=2.4):
        ncol = -(-ntor // nband)            # torsion columns per band
        rows = -(-nrow // nblock)           # population rows per block
        grid_w = nblock * (1.3 + panel_w * ncol) + (nblock - 1) * gap
        grid_h = nband * panel_h * rows + (nband - 1) * gap
        mol_w = min(max(4.5, 0.2 * grid_w), max(3.0, 0.9 * grid_h))   # larger molecules get more room
        return ncol, rows, grid_w, grid_h, mol_w

    def fit(nband, nblock):
        """Panel height (1.3-2.6 in) and, if the tallest panels are still too wide a figure, a
        narrower panel width (2.4 down to 1.6 in), that bring this split closest to 16:9."""
        panel_h, panel_w = 1.55, 2.4
        for _ in range(6):   # mol_w depends on the grid size, so settle in a few passes
            ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
            want_h = (mol_w + grid_w) / target - 1.0
            panel_h = float(np.clip((want_h - (nband - 1) * gap) / (nband * rows), 1.3, 2.6))
            ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
            if (mol_w + grid_w) / (grid_h + 1.0) > target:   # still too wide: narrow the panels
                want_grid_w = target * (grid_h + 1.0) - mol_w
                per_block = (want_grid_w - (nblock - 1) * gap) / nblock - 1.3
                panel_w = float(np.clip(per_block / ncol, 1.6, 2.4))
        ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
        return panel_h, panel_w, (mol_w + grid_w) / (grid_h + 1.0)

    # one split direction at a time; ties go to the simpler layout (fewer splits come first)
    splits = [(1, 1)] + [(b, 1) for b in (2, 3) if b <= ntor] + [(1, k) for k in (2, 3) if k <= nrow]
    nband, nblock = min(splits, key=lambda bk: round(abs(np.log(fit(*bk)[2] / target)), 3))
    panel_h, panel_w, _ = fit(nband, nblock)
    ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
    # text scales with the figure: unchanged up to ~20 in wide, larger on big grids so it stays
    # readable when the figure is shrunk to page width
    fs = max(1.0, (mol_w + grid_w) / 20.0)
    fig = plt.figure(figsize=(mol_w + grid_w, grid_h + 1.0), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=[mol_w, grid_w])
    mol_ax = fig.add_subplot(outer[0])
    png = molecule_png(mol, quads, colours, size=(int(220 * mol_w), int(220 * mol_w)))
    mol_ax.imshow(mpimg.imread(io.BytesIO(png), format="png"))
    mol_ax.axis("off")
    mol_ax.set_title("t<k> marks the central bond of torsion k", fontsize=10 * fs)
    cells = outer[1].subgridspec(nband, nblock, hspace=gap / (panel_h * rows),
                                 wspace=gap / (1.3 + panel_w * ncol))
    grids = [[cells[b, k].subgridspec(rows, ncol) for k in range(nblock)] for b in range(nband)]
    best = best or {}
    labels = [f"{top_label}\n(noise {100 * noise:.1f}%)"] + (list(row_labels) if row_labels else [
        f"group {g}\n({100 * float(masses[i]):.1f}%)"
        + "".join(f"\n◆ crystal {c}" for c in best.get(g, [])) for i, g in enumerate(groups)])

    edges = np.linspace(-np.pi, np.pi, 49)
    for t, (name, colour) in enumerate(zip(names, colours)):
        band, col = divmod(t, ncol)
        for pop in range(nrow):
            block, row = divmod(pop, rows)
            ax = fig.add_subplot(grids[band][block][row, col])
            values = tors[:, t] if pop == 0 else aligned[groups[pop - 1]][:, t]
            if log_floor:   # fraction of the row's frames per bin, log axis cut at log_floor
                ax.hist(values, bins=edges, weights=np.full(len(values), 1.0 / max(len(values), 1)),
                        color=colour, alpha=0.9 if pop == 0 else 0.75, histtype="stepfilled",
                        edgecolor=colour, bottom=log_floor)
                ax.set_yscale("log")
                ax.set_ylim(log_floor, 1.0)
            else:
                ax.hist(values, bins=edges, density=True, color=colour,
                        alpha=0.9 if pop == 0 else 0.75,
                        histtype="stepfilled", edgecolor=colour)
            if pop == 0:
                row_refs = refs if top_refs is None else top_refs
            else:
                row_refs = cluster_refs[groups[pop - 1]] if cluster_refs else refs
            draw_reference_lines(ax, row_refs, t, with_label=(pop == 0 and t == 0))
            ax.set_xlim(-np.pi, np.pi)
            ax.set_xticks([-np.pi, 0, np.pi], ["-π", "0", "π"])
            if log_floor and col == 0:
                ax.set_yticks([log_floor, 1e-3, 1e-1], [f"{log_floor:.0e}", "1e-3", "1e-1"])
                ax.minorticks_off()
            else:
                ax.set_yticks([])
                ax.minorticks_off()
            ax.tick_params(labelsize=7 * fs)
            if row == 0:
                ax.set_title(f"t{t}\n{name}", color=colour, fontsize=8 * fs)
            if pop == 0:
                ax.set_facecolor("#f4f4f4")
            elif best.get(groups[pop - 1]):
                for spine in ax.spines.values():
                    spine.set_edgecolor(BEST_MATCH_COLOUR)
                    spine.set_linewidth(3.0)
            if col == 0:
                ax.set_ylabel(labels[pop], fontsize=8 * fs, rotation=0, ha="right", va="center",
                              color="#8a6000" if pop and best.get(groups[pop - 1]) else "black")
            if pop == 0 and t == 0 and refs:
                ax.legend(fontsize=6 * fs, loc="upper left", frameon=False)
    fig.suptitle(title, fontsize=12 * fs)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def build_fit(tors, cmd, mol, torsions, *, allow_single_cluster=False, max_frames=None,
              max_matches=None, mass_floor=None, k=None, min_vote_fraction=None):
    """Shu-Yu's symmetry-first t_hdbscan (MD-tools main, docs/tutorial/paracetamol/analysis/
    symmetry-alignment.ipynb): TorsionDefinitions on the build SDF, the atom mapping checked
    against the trajectory topology, and the trajectory itself so the symmetry images of the
    selected torsions can be recomputed. Returns (fit, mapping, why)."""
    definitions = TorsionDefinitions(quads=tuple(torsions.values()), names=tuple(torsions),
                                     indexing="zero-based", units="radian")
    bonds = [(b[0].index, b[1].index) for b in cmd.topology.bonds]
    if not bonds:
        sys.exit("the topology has no bonds - use one with CONECT records / bonds")
    mapping, why = verify_atom_mapping(mol, bonds, [a.element.symbol for a in cmd.topology.atoms])
    extra = {} if max_frames is None else {"max_precomputed_frames": int(max_frames)}
    if max_matches is not None:
        extra["max_matches"] = int(max_matches)
    if mass_floor is not None:          # HDBSCAN min_cluster_size as a fraction of the frames
        extra["mass_floor"] = float(mass_floor)
    if k is not None:                   # voters: nearest clustered training frames per frame
        extra["k"] = int(k)
    if min_vote_fraction is not None:   # share of the k voters that must agree on the label
        extra["min_vote_fraction"] = float(min_vote_fraction)
    fit = t_hdbscan(tors, molecule=mol, torsion_definitions=definitions, coordinates=cmd,
                    atom_mapping=mapping, allow_single_cluster=allow_single_cluster, **extra).fit()
    return fit, mapping, why


def gauge_of(args, names):
    """The `gauge=` for align_all from --gauge (a torsion name or 0-based position) and
    --gauge-width (radians, a half-turn by default, centred on the reference), or None."""
    if not args.gauge:
        return None
    column = names[int(args.gauge)] if args.gauge.isdigit() else args.gauge
    if column not in names:
        sys.exit(f"--gauge {args.gauge!r}: not one of the torsions {names}")
    return (column, float(args.gauge_width))


def gauge_crystal(rows, reference_row, al, desc, fallback):
    """Crystal line in the gauge's convention (drawing only): of the crystal's atom assignments
    (`rows`, from reference_assignments), those that put the gauge column inside al's window, and
    of these the one nearest `reference_row`. `fallback` (the labelling chosen without a gauge) is
    kept when there is no gauge or no assignment lands in the window (a crystal on the cut)."""
    if "gauge_window" not in al:
        return fallback
    col, lo, hi = al["gauge_window"]
    offset = np.mod(rows[:, col] - lo, 2.0 * np.pi)
    inside = rows[np.isfinite(rows[:, col]) & (offset < hi - lo)]
    if not len(inside):
        return fallback
    return desc.nearest_labelling(inside, reference_row)["angles"]


def report_gauge(al, fit, tors_col, label=""):
    """Shu-Yu's checks: frames outside the window (on the cut), frames near it, flips used."""
    ok = al["gauge_satisfied"]
    near = al["gauge_edge_distance"] < 0.1
    col, lo, hi = al["gauge_window"]
    flipped = np.abs(np.angle(np.exp(1j * (al["angles"][:, col] - fit.closed_theta_[:, col])))) > 1.0
    print(f"{label}gauge window [{np.rad2deg(lo):.0f}, {np.rad2deg(hi):.0f}) deg on column {col}: "
          f"{100 * ok.mean():.2f}% in the window ({int((~ok).sum())} frames on the cut), "
          f"{100 * near.mean():.2f}% within 0.1 rad of a cut, flip applied to {100 * flipped.mean():.1f}%")
    for c in range(fit.n_clusters_):
        m = fit.labels_ == c
        print(f"{label}  cluster {group_name(c)}: flip applied {100 * flipped[m].mean():.1f}%, "
              f"not in window {int((~ok[m]).sum())}, within 0.1 rad of a cut {100 * near[m].mean():.1f}%")


def recursive_split(top_fit, tors, cmd, mol, torsions, X, operations, args):
    """Re-cluster every cluster whose internal (within-cluster, symmetry-corrected heavy-atom) RMSD
    is above args.recursive_max_rmsd on its own frames, until none is or a cluster cannot be split
    (no cluster / one sub-cluster) or args.recursive_max_depth is reached. Same fit as the top level
    (build_fit, same vote); the mass floor is a fraction of the frames re-clustered, raised to at
    least k / n so the vote cannot empty a cluster. Frames a re-clustering leaves unassigned are
    unassigned in the result. Names follow the tree: A, A.A, A.B, ...

    Returns {"nodes": name -> {fit, al, idx, children}, "leaves": [(leaf, node, local cluster)]
    sorted by size, "labels": final leaf index per frame (-1 unassigned), "tree": rows}."""
    cap = args.rmsd_max_per_cluster or None

    def internal(fit, idx):
        r = intercluster_rmsd(X[idx], fit.labels_, operations, max_per_cluster=cap, workers=args.rmsd_workers)
        return {c: float(r["mean"][i, i]) for i, c in enumerate(r["clusters"])}

    def refit(idx):
        floor = max(args.mass_floor or 0.01, (args.k or 20) / len(idx))
        for _ in range(4):
            try:
                sub, _, _ = build_fit(tors[idx], cmd[idx], mol, torsions, max_frames=args.max_frames,
                                      max_matches=args.max_matches, mass_floor=floor, k=args.k,
                                      min_vote_fraction=args.min_vote_fraction)
                return sub, sub.align_all(per_cluster=True, gauge=gauge_of(args, list(torsions))), floor
            except RuntimeError as e:
                if "no draw produced a cluster" in str(e):
                    return None, None, floor
                raise
            except ValueError as e:            # the vote emptied a cluster: raise the floor
                if "has no members" not in str(e):
                    raise
                floor *= 2
        return None, None, floor

    n = len(tors)
    nodes, leaves, tree = {}, [], []
    queue = [("", top_fit, top_fit.align_all(per_cluster=True, gauge=gauge_of(args, list(torsions))),
              np.arange(n), 0)]
    while queue:
        name, fit, al, idx, depth = queue.pop(0)
        within = internal(fit, idx)
        nodes[name] = {"fit": fit, "al": al, "idx": idx, "children": {}}
        for c in range(fit.n_clusters_):
            child = (name + "." if name else "") + group_name(c)
            members = idx[fit.labels_ == c]
            nodes[name]["children"][c] = child
            row = {"cluster": child, "parent": name, "depth": depth + 1, "frames": len(members),
                   "mass_of_ensemble": round(len(members) / n, 5), "internal_rmsd_A": round(within[c], 3)}
            outcome = "leaf"
            if within[c] > args.recursive_max_rmsd:
                if depth + 1 >= args.recursive_max_depth:
                    outcome = "leaf (max depth)"
                else:
                    sub, sub_al, floor = refit(members)
                    if sub is not None and sub.n_clusters_ >= 2:
                        queue.append((child, sub, sub_al, members, depth + 1))
                        outcome = f"split (mass floor {floor:.4f})"
                    else:
                        outcome = "leaf (cannot split: " + ("one state)" if sub is None else "one sub-cluster)")
            row["outcome"] = outcome
            tree.append(row)
            print(f"{'  ' * depth}{child}: {100 * len(members) / n:5.1f}% of the frames, internal RMSD "
                  f"{within[c]:.2f} A -> {outcome}")
            if not outcome.startswith("split"):
                leaves.append((child, name, c))
    size = lambda leaf: int(np.sum(nodes[leaf[1]]["fit"].labels_ == leaf[2]))
    leaves.sort(key=lambda leaf: -size(leaf))
    labels = np.full(n, -1)
    for i, (leaf, node, c) in enumerate(leaves):
        nd = nodes[node]
        labels[nd["idx"][nd["fit"].labels_ == c]] = i
    return {"nodes": nodes, "leaves": leaves, "labels": labels, "tree": tree}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True, help="solute trajectory of the physical state (state 0)")
    ap.add_argument("--top", required=True, help="matching topology WITH bonds (e.g. built.solute.pdb)")
    ap.add_argument("--sdf", required=True, help="ligand SDF from the build (bond orders)")
    ap.add_argument("--torsion", action="append", default=None,
                    help="name:i,j,k,l (repeat); default = tutorial paracetamol set")
    ap.add_argument("--torsions-file", action="append", default=None,
                    help="torsions.txt from select_torsions.py (repeatable; lines 'i j k l ...  # names')")
    ap.add_argument("--use-torsions", default=None,
                    help="comma-separated positions in the final torsion list to keep")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--frames-file", default=None,
                    help="re-cluster a SUBSET: text file of frame indices (one per line, the 'frame' "
                         "column of a previous run's labels.csv, i.e. indices after --stride), e.g. the "
                         "frames of one group; use the same --traj and --stride as that run")
    ap.add_argument("--allow-single-cluster", action="store_true",
                    help="allow_single_cluster=True: opt-in for an ensemble that may be ONE state "
                         "(MD-tools' default, and this script's, is False since 0.6.5.1)")
    ap.add_argument("--no-single-cluster", action="store_true",
                    help="allow_single_cluster=False; now the default, kept so old commands still run")
    ap.add_argument("--max-matches", type=int, default=None,
                    help="cap on the symmetry (automorphism) enumeration (MD-tools default 200000); "
                         "molecules with several methyl/CF3 groups exceed it and the fit refuses an "
                         "incomplete search -- raise it, e.g. 8000000")
    ap.add_argument("--mass-floor", type=float, default=None,
                    help="smallest cluster as a fraction of the fitted frames (HDBSCAN min_cluster_size; "
                         "MD-tools default 0.01). min_samples and the vote are unchanged")
    ap.add_argument("--k", type=int, default=None,
                    help="number of voters per frame (MD-tools default 20); a cluster needs about k "
                         "training frames to win a vote, so a low --mass-floor needs a low --k")
    ap.add_argument("--min-vote-fraction", type=float, default=None,
                    help="share of the k voters that must agree (MD-tools default 0.90)")
    ap.add_argument("--max-frames", type=int, default=None,
                    help="raise t_hdbscan's max_precomputed_frames (default 10000): the symmetry-aware "
                         "distance is a dense N x N matrix, 8 N^2 bytes per copy, about three copies")
    ap.add_argument("--name", default=None, help="compound name for the figure titles")
    ap.add_argument("--rest2-config", default=None,
                    help="REST2.config of the run: number of replicas and tau_max for the title")
    ap.add_argument("--duration-ns", type=float, default=None,
                    help="simulated time for the title (default: from the trajectory time stamps)")
    ap.add_argument("--ref-sdf", action="append", default=None,
                    help="crystal (bioactive) ligand SDF; repeatable (e.g. one per chain)")
    ap.add_argument("--pair-2d", action="append", default=None, metavar="X,Y",
                    help="draw the joint distribution of torsion Y (vertical axis) against torsion X "
                         "(horizontal), 0-based positions; repeatable. Default: 3,4")
    ap.add_argument("--gauge", default=None,
                    help="figures: write this torsion (name or 0-based position) in a window --gauge-width "
                         "wide centred on the density peak, the symmetry partners following it (MD-tools "
                         "0.6.5.2 gauge, e.g. one ring torsion of a para-phenyl). Clustering unchanged")
    ap.add_argument("--gauge-width", type=float, default=np.pi, help="window width, radians (default pi)")
    ap.add_argument("--recursive-max-rmsd", type=float, default=None,
                    help="re-cluster every cluster whose internal symmetry-corrected heavy-atom RMSD is "
                         "above this (A, e.g. 2.0) on its own frames, recursively; the final clusters "
                         "replace the top-level ones in every output (tree.csv records the hierarchy)")
    ap.add_argument("--recursive-max-depth", type=int, default=4)
    ap.add_argument("--rmsd-max-per-cluster", type=int, default=1500,
                    help="frames per cluster for the inter-cluster RMSD (Shu-Yu's c(RGDfV) setting, "
                         "1500); 0 = every frame (the paracetamol tutorial)")
    ap.add_argument("--rmsd-workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--out", default="tutorial-clustering")
    args = ap.parse_args()
    if args.allow_single_cluster and args.no_single_cluster:
        sys.exit("--allow-single-cluster and --no-single-cluster contradict each other")
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    P = lambda name: os.path.join(out, name)

    torsions = load_torsions(args.torsions_file, args.torsion, args.use_torsions)
    names, quads = list(torsions), list(torsions.values())
    label = f"{args.name}: " if args.name else ""
    print("md-tools", importlib.metadata.version("md-tools"), "from",
          os.path.dirname(md_tools.__file__))

    # ------------------------------------------------- 1. symmetry-first clustering
    print("== 1. Symmetry-first clustering (Shu-Yu, MD-tools main) ==")
    cmd = md.load(args.traj, top=args.top, stride=args.stride)
    frame_ids = np.arange(cmd.n_frames)
    if args.frames_file:
        frame_ids = np.loadtxt(args.frames_file, dtype=int, ndmin=1)
        if frame_ids.size == 0 or frame_ids.max() >= cmd.n_frames or frame_ids.min() < 0:
            sys.exit(f"{args.frames_file}: indices must lie in 0..{cmd.n_frames - 1} (same --traj/--stride)")
        cmd = cmd[frame_ids]
        print(f"subset: {cmd.n_frames} frames from {args.frames_file}")
    for name, q in torsions.items():
        print(f"  {name:8s} {q}: " + " - ".join(str(cmd.topology.atom(i)) for i in q))
    tors = md.compute_dihedrals(cmd, quads)
    print(f"{tors.shape[0]} frames, {tors.shape[1]} torsions")
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    if mol is None:
        sys.exit(f"RDKit could not read {args.sdf}")
    with open(P("torsion_cv.json"), "w") as f:
        json.dump({"schema_version": 1, "indexing": "zero-based", "units": "radian",
                   "torsions": [{"name": n, "atoms": list(q)} for n, q in torsions.items()]},
                  f, indent=2)

    asc = args.allow_single_cluster
    try:
        fit, mapping, why = build_fit(tors, cmd, mol, torsions, allow_single_cluster=asc,
                                      max_frames=args.max_frames, max_matches=args.max_matches,
                                      mass_floor=args.mass_floor, k=args.k,
                                      min_vote_fraction=args.min_vote_fraction)
    except RuntimeError as e:
        if "no draw produced a cluster" not in str(e) or asc:
            raise
        # allow_single_cluster=False and nothing admissible: EOM refused the root, which is what a
        # ONE-state ensemble gives (MD-tools 0.6.5.1 release note). Say so instead of a traceback.
        msg = ("NO CLUSTER with allow_single_cluster=False: no basin above the 1 % mass floor once "
               "the root of the condensed tree is refused -- the ensemble is most likely ONE state. "
               "Rerun with --allow-single-cluster. MD-tools: " + str(e))
        print(msg)
        with open(P("NO_CLUSTER.txt"), "w") as f:
            f.write(msg + "\n")
        sys.exit(3)
    print("atom mapping:", "identity" if mapping == list(range(mol.GetNumAtoms())) else mapping)
    print(" ", why)
    summary = fit.summary()
    desc = fit.descriptor_
    sym = summary["symmetry"]
    print(f"allow_single_cluster = {asc}; {summary['n_fit_frames']} frames fitted, "
          f"min_cluster_size {summary['settings']['min_cluster_size']}, "
          f"min_samples {summary['settings']['min_samples']}")
    print(f"symmetry operations acting on the torsions: {len(desc.operations)}")
    print("descriptor columns (selected + symmetry images):")
    for n, q, w, o in zip(desc.names, desc.quads, desc.column_weights, desc.origin):
        origin = "selected" if o["kind"] == "selected" else f"image of {o['image_of']}"
        print(f"  {n:28s} {'-'.join(map(str, q)):14s} weight {w:.3f}  {origin}")
    masses = {i: float(fit.cluster_mass_[i]) for i in range(fit.n_clusters_)}
    noise = float(fit.unassigned_mass_)
    print(f"n_clusters = {fit.n_clusters_}")
    # per-frame cluster labels, so a group can be re-clustered on its own (--frames-file)
    with open(P("labels.csv"), "w") as f:
        f.write("frame,time_ps,label\n")
        for fid, t_ps, lab in zip(frame_ids, cmd.time, fit.labels_):
            name = group_name(lab) if lab >= 0 else ("noise" if lab == -2 else "ambiguous")
            f.write(f"{int(fid)},{float(t_ps):.1f},{name}\n")
    print(f"masses     = { {group_name(k): round(v, 4) for k, v in masses.items()} }")
    print(f"density noise (-2) = {fit.noise_mass_:.4f}, ambiguous vote (-1) = "
          f"{fit.ambiguous_mass_:.4f}, unassigned total = {noise:.4f}")
    print("Caveats:")
    for caveat in summary["caveats"]:
        print("  -", caveat)
    with open(P("summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # heavy-atom coordinates and symmetry operations for the RMSDs (enumerated once: slow for 4j7d)
    heavy = heavy_atoms(mol)
    sym_kw = {} if args.max_matches is None else {"max_matches": int(args.max_matches)}
    operations = operations_on_atoms(enumerate_symmetry(mol, **sym_kw), heavy)
    X = cmd.xyz[:, np.asarray(mapping)[heavy], :] * 10.0      # nm -> A, heavy atoms in the molecule's order
    cap = args.rmsd_max_per_cluster or None

    rec = None
    if args.recursive_max_rmsd and fit.n_clusters_:
        print(f"\n== Recursive splitting: re-cluster while internal RMSD > {args.recursive_max_rmsd} A ==")
        rec = recursive_split(fit, tors, cmd, mol, torsions, X, operations, args)
        pd.DataFrame(rec["tree"]).to_csv(P("tree.csv"), index=False)
        with open(P("labels.csv"), "w") as f:          # final (leaf) cluster per frame
            f.write("frame,time_ps,label\n")
            for fid, t_ps, lab in zip(frame_ids, cmd.time, rec["labels"]):
                f.write(f"{int(fid)},{float(t_ps):.1f},{rec['leaves'][lab][0] if lab >= 0 else 'unassigned'}\n")
        print("final clusters:", ", ".join(f"{leaf} {100 * np.mean(rec['labels'] == i):.1f}%"
                                         for i, (leaf, _, _) in enumerate(rec["leaves"])),
              f"| unassigned {100 * np.mean(rec['labels'] < 0):.1f}%")

    # crystal ligands: every atom assignment, on all descriptor columns (Shu-Yu's tutorial)
    refs, crystal_rows = {}, {}
    for path in (args.ref_sdf or []):
        ref_label = os.path.splitext(os.path.basename(path))[0]
        crystal_rows[ref_label] = reference_assignments(path, mol, desc.quads)
        refs[ref_label] = crystal_rows_to_refs(crystal_rows[ref_label], len(names))
        print(f"crystal reference {ref_label}: {len(crystal_rows[ref_label])} assignment(s)")
        for name, values in zip(names, refs[ref_label]):
            shown = ", ".join(f"{np.rad2deg(v):.1f}" for v in values) or "n/a (through H)"
            print(f"  {name:20s} {shown} deg")

    # ------------------------- 2. one labelling for the figures (two-level Dirichlet alignment)
    print("\n== 2. Alignment: align_all(per_cluster=True) ==")
    gauge = gauge_of(args, names)
    al = fit.align_all(per_cluster=True, gauge=gauge)
    if gauge:
        report_gauge(al, fit, names.index(gauge[0]))
    groups = {group_name(c): c for c in range(fit.n_clusters_)}
    aligned ={g: al["selected"][fit.labels_ == c] for g, c in groups.items()}
    raw = {g: tors[fit.labels_ == c] for g, c in groups.items()}
    print(f"reference: frame {al['reference']} (density peak of cluster A)")
    for g, c in groups.items():
        m = fit.labels_ == c
        print(f"  cluster {g}: density peak frame {int(al['cluster_references'][c])}, "
              f"relabelled {100 * float((al['operations'][m] != 0).mean()):.1f}%, "
              f"margin < 0.1 for {100 * float((al['margin'][m] < 0.1).mean()):.1f}%")
    reps = fit.representatives()
    chosen = {g: int(reps[c]) for g, c in groups.items()}
    peaks = {g: int(al["cluster_references"][c]) for g, c in groups.items()}
    if rec:                       # the final (leaf) clusters replace the top-level ones from here on
        groups, aligned, raw, chosen, peaks, node_reps = {}, {}, {}, {}, {}, {}
        for i, (leaf, node, c) in enumerate(rec["leaves"]):
            nd = rec["nodes"][node]
            m = nd["fit"].labels_ == c
            if node not in node_reps:
                node_reps[node] = nd["fit"].representatives()
            groups[leaf] = i
            aligned[leaf] = nd["al"]["selected"][m]
            raw[leaf] = tors[nd["idx"][m]]
            chosen[leaf] = int(nd["idx"][int(node_reps[node][c])])
            peaks[leaf] = int(nd["idx"][int(nd["al"]["cluster_references"][c])])
        masses = {i: float(np.mean(rec["labels"] == i)) for i in range(len(groups))}
        noise = float(np.mean(rec["labels"] < 0))

    # crystal -> cluster: per cluster the assignment nearest its density peak (fit.assign_reference)
    cluster_refs = {g: {} for g in groups}
    best_groups, crystal_notes = {}, []
    assignment_rows, assignment_json = [], {}
    for ref_label, rows in crystal_rows.items() if rec else ():
        # follow the tree: assign at the top level, then inside the chosen cluster's split, ...
        node, outs = "", {}
        for nname, nd in rec["nodes"].items():
            outs[nname] = nd["fit"].assign_reference(rows)
        while True:
            o = outs[node]
            b = o["best_cluster"]
            if b is None:
                final_best = None
                consistent = [rec["nodes"][node]["children"][c] for c in o["consistent_with"]]
                break
            child = rec["nodes"][node]["children"][b]
            if child in rec["nodes"]:
                node = child
                continue
            final_best, consistent = child, [child]
            break
        print(f"crystal {ref_label}: best cluster {final_best or 'none'}; cannot be told from {consistent}")
        if final_best:
            best_groups.setdefault(final_best, []).append(ref_label)
        else:
            crystal_notes.append(f"crystal {ref_label}: no single best group (cannot tell "
                                 f"{', '.join(consistent)} apart)")
        assignment_json[ref_label] = {"best_cluster": final_best, "consistent_with": consistent,
                                      "assigned_through": node or "(top level)"}
        for leaf, nname, c in rec["leaves"]:
            v = outs[nname]["clusters"][c]
            nal = rec["nodes"][nname]["al"]
            line = gauge_crystal(rows, nal["angles"][int(nal["cluster_references"][c])], nal, desc, v["angles"])
            cluster_refs[leaf][ref_label] = [[float(a)] if np.isfinite(a) else [] for a in line[:len(names)]]
            assignment_rows.append({"crystal": ref_label, "cluster": leaf, "assignment chosen": v["assignment"],
                                    "distance to the density peak": round(v["distance_to_peak"], 3),
                                    "nearest frame": round(v["nearest_member_distance"], 3),
                                    "frames nearer the peak than the crystal": round(v["peak_percentile"], 3),
                                    "best cluster": final_best or "none"})
    for ref_label, rows in crystal_rows.items() if not rec else ():
        out_ref = fit.assign_reference(rows)
        best = out_ref["best_cluster"]
        print(f"crystal {ref_label}: columns compared {out_ref['columns_used']}; best cluster "
              f"{group_name(best) if best is not None else 'none'}; cannot be told from "
              f"{[group_name(c) for c in out_ref['consistent_with']]}")
        if best is not None:
            best_groups.setdefault(group_name(best), []).append(ref_label)
        else:
            crystal_notes.append(f"crystal {ref_label}: no single best group (cannot tell "
                                 f"{', '.join(group_name(c) for c in out_ref['consistent_with'])} apart)")
        for pair in out_ref["indistinguishable"]:
            print(f"   clusters {[group_name(c) for c in pair['clusters']]}: density peaks "
                  f"{pair['peak_gap']:.3f} apart on these columns, inside a cluster's own spread "
                  f"({pair['spread']:.3f})")
        for g, c in groups.items():
            v = out_ref["clusters"][c]
            line = gauge_crystal(rows, al["angles"][int(al["cluster_references"][c])], al, desc, v["angles"])
            cluster_refs[g][ref_label] = [[float(a)] if np.isfinite(a) else [] for a in line[:len(names)]]
            assignment_rows.append({"crystal": ref_label, "cluster": g, "assignment chosen": v["assignment"],
                                    "distance to the density peak": round(v["distance_to_peak"], 3),
                                    "nearest frame": round(v["nearest_member_distance"], 3),
                                    "frames nearer the peak than the crystal": round(v["peak_percentile"], 3),
                                    "best cluster": group_name(best) if best is not None else "none"})
        assignment_json[ref_label] = {
            "best_cluster": group_name(best) if best is not None else None,
            "consistent_with": [group_name(c) for c in out_ref["consistent_with"]],
            "indistinguishable": [{"clusters": [group_name(c) for c in p["clusters"]],
                                   "peak_gap": p["peak_gap"], "spread": p["spread"]}
                                  for p in out_ref["indistinguishable"]],
            "columns_used": out_ref["columns_used"],
            "clusters": {group_name(c): {k: (np.rad2deg(v[k]).round(2).tolist() if k == "angles" else v[k])
                                         for k in v} for c, v in out_ref["clusters"].items()}}
    if assignment_rows:
        pd.DataFrame(assignment_rows).to_csv(P("crystal_assignment.csv"), index=False)
        with open(P("crystal_assignment.json"), "w") as f:
            json.dump(assignment_json, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))

    # ---------------------------------------------------------------- 3. summary + figures
    print("\n== 3. Summary ==")
    for i, (g, c) in enumerate(groups.items()):
        print(f"group {g}: cluster {c}  mass {masses[i]:.4f}  representative frame {chosen[g]}")
    group_list = list(groups)
    with open(P("merged_groups.json"), "w") as f:     # same file name as before, for downstream
        json.dump({"route": summary["route"], "groups": {g: [c] for g, c in groups.items()},
                   "masses": {g: masses[i] for i, g in enumerate(group_list)},
                   "noise_mass": float(fit.noise_mass_), "ambiguous_mass": float(fit.ambiguous_mass_),
                   "representatives": chosen, "alignment_reference": int(al["reference"]),
                   "cluster_density_peaks": peaks}, f, indent=2)

    # (group_torsions.png, all groups overlaid per torsion, dropped 2026-10-26: not informative)

    n_rep, tau_max = read_rest2_config(args.rest2_config)
    duration = args.duration_ns if args.duration_ns is not None else simulated_ns(cmd, args.stride)
    run_info = [f"{duration:g} ns" if duration else None,
                f"{n_rep} replicas" if n_rep else None,
                f"τ_max {tau_max}" if tau_max is not None else None,
                f"stride {args.stride}, {tors.shape[0]} frames"]
    head = f"{args.name or 'ligand'} — REST2 " + ", ".join(x for x in run_info if x)
    crystal = "; black lines = crystal (bioactive) conformation" if refs else ""
    marker = (f"; gold frame = group the crystal is assigned to" if best_groups else "")
    notes = ("\n" + "; ".join(crystal_notes)) if crystal_notes else ""
    # top row in ONE labelling as well (every frame relabelled by align_all, the same alignment as
    # the cluster rows), with each crystal in the labelling nearest the alignment reference -- so a
    # folded torsion (e.g. a carboxylate) shows one peak; the as-stored view is the _unaligned figure
    top_refs = {}
    for ref_label, rows in crystal_rows.items():
        near = desc.nearest_labelling(rows, al["angles"][int(al["reference"])])
        near = {"angles": gauge_crystal(rows, al["angles"][int(al["reference"])], al, desc, near["angles"])}
        top_refs[ref_label] = [[float(a)] if np.isfinite(a) else [] for a in near["angles"][:len(names)]]
    unassigned = (f"unassigned {100 * noise:.1f}% (density noise {100 * fit.noise_mass_:.1f}%, "
                  f"ambiguous {100 * fit.ambiguous_mass_:.1f}%)")
    if gauge:
        crystal += f"; gauge: {gauge[0]} in a half-turn window on the density peak"
    if rec:
        unassigned = (f"recursive: split while internal RMSD > {args.recursive_max_rmsd} Å, "
                      f"unassigned {100 * noise:.1f}% (top level {100 * float(fit.unassigned_mass_):.1f}%)")
    draw_overview(P("clusters_overview.png"), mol=mol, names=names, quads=quads, tors=al["selected"],
                  aligned=aligned, masses=masses, noise=noise, refs=refs, cluster_refs=cluster_refs,
                  top_refs=top_refs, top_label="all frames,\none labelling", best=best_groups,
                  title=f"{head}\n{len(groups)} symmetry-first clusters, {unassigned}{crystal}{marker}{notes}")

    # The same overview with each cluster's torsions AS STORED (no relabelling): a cluster that
    # holds one conformer in two labellings shows both images here (e.g. peaks shifted by pi).
    draw_overview(P("clusters_overview_unaligned.png"), mol=mol, names=names, quads=quads,
                  tors=tors, aligned=raw, masses=masses, noise=noise, refs=refs,
                  top_label="all frames,\nas stored", best=best_groups,
                  title=f"{head}\n{len(groups)} symmetry-first clusters, torsions AS STORED "
                        f"(not relabelled), {unassigned}{crystal}")

    # ------------------------------------------------- 2-D torsion marginals
    pairs = [tuple(int(v) for v in p.split(",")) for p in (args.pair_2d or ["3,4"])]
    for ix, iy in pairs:               # 0-based torsion positions: x axis, y axis
        if len(names) <= max(ix, iy):
            print(f"  (skipping 2-D plot {ix},{iy}: only {len(names)} torsions)")
            continue
        panels = 1 + len(groups)
        ncol2 = min(3, panels)
        nrow2 = int(np.ceil(panels / ncol2))
        fig2, axes2 = plt.subplots(nrow2, ncol2, figsize=(4.5 * ncol2, 4.1 * nrow2),
                                   layout="constrained", squeeze=False)
        for extra in axes2.flat[panels:]:
            extra.axis("off")
        lim = (-np.pi, np.pi, -np.pi, np.pi)

        def _hex2d(axis, xs, ys, title):
            hb = axis.hexbin(xs, ys, gridsize=55, mincnt=1, cmap="viridis",
                             bins="log", extent=lim)
            fig2.colorbar(hb, ax=axis, label="log10(count)")
            axis.set_title(title)
            axis.set_xlabel(f"{names[ix]} (rad)")
            axis.set_ylabel(f"{names[iy]} (rad)")

        _hex2d(axes2.flat[0], tors[:, ix], tors[:, iy], "all frames")
        for k, (g, values) in enumerate(aligned.items()):
            _hex2d(axes2.flat[k + 1], values[:, ix], values[:, iy],
                   f"group {g} ({100 * masses[k]:.1f}%)")
        fig2.suptitle(f"{label}joint distribution: {names[iy]} as a function of {names[ix]}", fontsize=11)
        fig2.savefig(P(f"torsions_{ix}_{iy}_2d.png"), dpi=150)
        plt.close(fig2)

    rows = []
    for i, (g, c) in enumerate(groups.items()):
        values = aligned[g]
        row = {"cluster_id": g, "original_clusters": str(c)}
        for column, torsion in enumerate(names):
            row[f"{torsion}_avg"] = np.rad2deg(circular_mean(values[:, column]))
            row[f"{torsion}_std"] = np.rad2deg(circular_std(values[:, column]))
        row["n_samples"] = int(values.shape[0])
        row["mass"] = masses[i]
        rows.append(row)
    df = pd.DataFrame(rows).round(3)
    print(df.to_string(index=False))
    df.to_csv(P("groups.csv"), index=False)

    # ------------------------------------------- inter-cluster RMSD, symmetry-corrected (0.6.5.1)
    if groups:
        final_labels = rec["labels"] if rec else fit.labels_
        cname = (lambda c: rec["leaves"][c][0]) if rec else group_name
        corrected = intercluster_rmsd(X, final_labels, operations, max_per_cluster=cap, workers=args.rmsd_workers)
        by_index = intercluster_rmsd(X, final_labels, operations[:1], max_per_cluster=cap, workers=args.rmsd_workers)
        names_c = [cname(c) for c in corrected["clusters"]]

        def as_table(r):
            return pd.DataFrame([[f"{r['mean'][i, j]:.2f} +- {r['std'][i, j]:.2f}" for j in range(len(names_c))]
                                 for i in range(len(names_c))], index=names_c, columns=names_c)
        print(f"\n{len(heavy)} heavy atoms, {corrected['n_operations']} symmetry operations, frames per "
              f"cluster {corrected['n_frames']} of {corrected['n_frames_in_cluster']}")
        print("symmetry-corrected RMSD (A), mean +- std:")
        print(as_table(corrected).to_string())
        print("by atom index (no symmetry), mean +- std:")
        print(as_table(by_index).to_string())
        rmsd_rows = [{"cluster_a": names_c[i], "cluster_b": names_c[j],
                      "rmsd_mean": corrected["mean"][i, j], "rmsd_std": corrected["std"][i, j],
                      "rmsd_by_index_mean": by_index["mean"][i, j], "rmsd_by_index_std": by_index["std"][i, j]}
                     for i in range(len(names_c)) for j in range(len(names_c))]
        pd.DataFrame(rmsd_rows).round(4).to_csv(P("intercluster_rmsd.csv"), index=False)
        with open(P("intercluster_rmsd.json"), "w") as f:
            json.dump({"heavy_atoms": len(heavy), "n_operations": corrected["n_operations"],
                       "max_per_cluster": cap,
                       "n_frames": {cname(c): n for c, n in corrected["n_frames"].items()},
                       "n_frames_in_cluster": {cname(c): n for c, n in corrected["n_frames_in_cluster"].items()},
                       "clusters": names_c, "corrected_mean": corrected["mean"].tolist(),
                       "corrected_std": corrected["std"].tolist(), "by_index_mean": by_index["mean"].tolist(),
                       "by_index_std": by_index["std"].tolist()}, f, indent=2)
        k = len(names_c)
        fig, axes = plt.subplots(1, 2, figsize=(4.5 + 0.9 * k, 2.2 + 0.8 * k), layout="constrained")
        vmax = max(np.nanmax(corrected["mean"]), np.nanmax(by_index["mean"]))
        for ax, (title, r) in zip(axes, (("symmetry-corrected", corrected), ("by atom index", by_index))):
            im = ax.imshow(r["mean"], cmap="viridis", vmin=0, vmax=vmax)
            ax.set_xticks(range(k), names_c)
            ax.set_yticks(range(k), names_c)
            for i in range(k):
                for j in range(k):
                    ax.text(j, i, f"{r['mean'][i, j]:.2f}\n+- {r['std'][i, j]:.2f}", ha="center", va="center",
                            color="w" if r["mean"][i, j] < 0.6 * vmax else "k", fontsize=9)
            ax.set_title(f"heavy-atom RMSD (Å), {title}", fontsize=10)
        fig.colorbar(im, ax=axes, shrink=0.8, label="Å")
        fig.suptitle(f"{label}between clusters (off-diagonal) and within (diagonal)", fontsize=11)
        fig.savefig(P("intercluster_rmsd.png"), dpi=150)
        plt.close(fig)

    keep = strip_nonpolar_hydrogens(cmd)
    try:
        info = draw_representative_structures(cmd, chosen, output=P("representatives.png"), atom_indices=keep)
        print(f"{info['n_atoms_drawn']} atoms drawn, superposed: {info['aligned']}")
    except Exception as e:
        print(f"  (could not draw representatives: {e})")
    for g, frame in chosen.items():
        cmd[frame].save_pdb(P(f"representative_{g}.pdb"))
    print(f"\nAll outputs in {out}/")


if __name__ == "__main__":
    main()
