#!/usr/bin/env python
"""
Folate torsions, REST2 v4 (Sage 2.2.1 torsions) vs FOL-REST2-mace (bespokefit MACE-OMOL torsions):
one figure per system (water, 4LRH, 3QL3), one panel per torsion, all state-0 frames, both
distributions overlaid (v4 blue filled, MACE orange line), crystal torsion(s) as black lines.

Why the torsions have to be matched first: the two packages number folate differently (v4 starts
O1 C1 O2..., MACE starts N1 C1 N2...) and their select_torsions.py lists differ in order and in the
outer atoms chosen. So every torsion is defined ONCE on v4's atoms and translated to MACE indices
through a heavy-atom graph match (element + connectivity only; bond orders, charges and H ignored,
so the N1-H vs 3H-pterin tautomer difference does not matter). The torsion set is the union of both
runs' torsions.txt + torsions_rigid.txt, one per central bond, ordered pteridine -> pABA -> glutamate.

Equivalent end atoms (carboxylate O's, a phenyl ring flip) are folded: a torsion whose last (or
first) atom has n topologically equivalent siblings is shown modulo 360/n degrees, so it does not
matter which of the equivalent atoms the graph match picked.

Each panel reports the overlap coefficient sum(min(p_v4, p_mace)) (1 = identical distributions,
0 = disjoint) and, as a noise reference, the overlap between the first and second half of the v4
run. Densities are normalised, so the different lengths (water 100 ns, complexes 50 ns) don't bias it.

Inputs (state 0, folate only, same atom order as build/FOL.sdf):
  water    <run>/water/REST2-run1/solute_state0_prod1.nc
  complex  <run>/<PDB>/analysis-shuyu/FOL_state0_prod1.nc

Needs numpy, scipy, matplotlib, rdkit (bsc-thesis/venv has them; no MD-tools needed).
  python compare_v4_mace_torsions.py [--thesis DIR] [--crystal DIR] [--out DIR] [--bin DEG]
"""
import argparse
import csv
import io
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from scipy.io import netcdf_file

RDLogger.DisableLog("rdApp.*")

HOME = os.path.expanduser("~")
SYSTEMS = ("water", "4LRH", "3QL3")
RUNS = (("v4", "FOL-REST2-v4", "#2a78d6"), ("MACE", "FOL-REST2-mace", "#eb6834"))
# v4 atom names (build/built.solute.pdb), from FOL-REST2-v4/scripts/fragments.py
FRAGMENTS = {"pteridine": "C11 C12 C13 C14 C15 C16 N3 N4 N5 N6 N7 O4",
             "pABA": "C4 C5 C6 C7 C8 C9 N2 C10 C3 O3",
             "glutamate": "N1 C2 C1 O1 O2 C17 C18 C19 O5 O6"}


def bare_graph(sdf):
    """Heavy atoms only, every bond single, no aromaticity/charge: element + connectivity."""
    mol = Chem.RemoveHs(Chem.MolFromMolFile(sdf, removeHs=False, sanitize=False), sanitize=False)
    rw = Chem.RWMol(mol)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def heavy_index(sdf):
    """Heavy-atom position in the bare graph -> atom index in the full molecule (= trajectory)."""
    mol = Chem.MolFromMolFile(sdf, removeHs=False, sanitize=False)
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]


def graph_map(ref, other):
    """ref heavy position -> other heavy position (one of the symmetry-equivalent matches)."""
    match = other.GetSubstructMatch(ref)
    if len(match) != ref.GetNumAtoms() or ref.GetNumAtoms() != other.GetNumAtoms():
        raise SystemExit("heavy-atom graphs of the two folate packages do not match")
    return dict(enumerate(match))


def read_torsions(run_dir):
    quads = []
    for name in ("torsions.txt", "torsions_rigid.txt"):
        path = f"{run_dir}/{name}"
        if os.path.exists(path):
            for line in open(path):
                if line.strip() and not line.startswith("#"):
                    quads.append(tuple(int(x) for x in line.split("#")[0].split()))
    return quads


def heavy_end(graph, j, k):
    """Heavy neighbour of j other than k (highest canonical rank, deterministic)."""
    rank = list(Chem.CanonicalRankAtoms(graph, breakTies=True))
    nb = [n.GetIdx() for n in graph.GetAtomWithIdx(j).GetNeighbors() if n.GetIdx() != k]
    return max(nb, key=lambda n: rank[n]) if nb else None


def fold_period(graph, quad):
    """n if an end atom has n topologically equivalent siblings (COO-: 2, ring flip: 2), else 1."""
    rank = list(Chem.CanonicalRankAtoms(graph, breakTies=False))
    period = 1
    for j, k, l in ((quad[2], quad[1], quad[3]), (quad[1], quad[2], quad[0])):
        sib = [n.GetIdx() for n in graph.GetAtomWithIdx(j).GetNeighbors()
               if n.GetIdx() != k and rank[n.GetIdx()] == rank[l]]
        period = max(period, len(sib))
    return period


def dihedral(xyz, quad):
    """Dihedral (deg) per frame; xyz is (frames, atoms, 3)."""
    p0, p1, p2, p3 = (xyz[:, i].astype(float) for i in quad)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 /= np.linalg.norm(b1, axis=1, keepdims=True)
    v = b0 - (b0 * b1).sum(1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(1, keepdims=True) * b1
    x = (v * w).sum(1)
    y = (np.cross(b1, v) * w).sum(1)
    return np.degrees(np.arctan2(y, x))


def fold(deg, period):
    """Map onto [-180/period, 180/period)."""
    span = 360.0 / period
    return (np.asarray(deg) + span / 2) % span - span / 2


def overlap(a, b, edges):
    pa, _ = np.histogram(a, bins=edges)
    pb, _ = np.histogram(b, bins=edges)
    return float(np.minimum(pa / pa.sum(), pb / pb.sum()).sum())


def load_coords(nc):
    f = netcdf_file(nc, "r", mmap=False)
    xyz = np.array(f.variables["coordinates"][:], dtype=np.float32)
    f.close()
    return xyz


def atom_names(pdb):
    return [line[12:16].strip() for line in open(pdb) if line.startswith(("ATOM", "HETATM"))]


def build_torsions(thesis):
    """Union of both runs' torsions on v4 heavy positions, one per central bond, plus the maps."""
    v4_sdf = f"{thesis}/FOL-REST2-v4/water/build/FOL.sdf"
    mc_sdf = f"{thesis}/FOL-REST2-mace/water/build/FOL.sdf"
    g_v4, g_mc = bare_graph(v4_sdf), bare_graph(mc_sdf)
    hv4, hmc = heavy_index(v4_sdf), heavy_index(mc_sdf)
    v4_to_mc = graph_map(g_v4, g_mc)
    pos_v4 = {a: i for i, a in enumerate(hv4)}
    pos_mc = {a: i for i, a in enumerate(hmc)}
    mc_to_v4 = {m: v for v, m in v4_to_mc.items()}

    by_bond, source = {}, {}
    for label, run, quads, pos, back in (("v4", "FOL-REST2-v4", None, pos_v4, None),
                                         ("MACE", "FOL-REST2-mace", None, pos_mc, mc_to_v4)):
        for q in read_torsions(f"{thesis}/{run}/water/REST2-run1"):
            p = [pos.get(a) for a in q]
            if p[1] is None or p[2] is None:
                continue  # central bond to H: not a heavy-atom torsion
            if back:
                p = [None if x is None else back[x] for x in p]
            if p[0] is None:
                p[0] = heavy_end(g_v4, p[1], p[2])
            if p[3] is None:
                p[3] = heavy_end(g_v4, p[2], p[1])
            key = frozenset(p[1:3])
            by_bond.setdefault(key, tuple(p))  # v4's own definition wins when both have it
            source.setdefault(key, set()).add(label)

    # order along the molecule: distance from the pteridine end of the longest path
    # (pteridine = the 6-rings containing >= 2 N; the pABA ring has none)
    dm = Chem.GetDistanceMatrix(g_v4)
    pter = {x for r in g_v4.GetRingInfo().AtomRings() if len(r) == 6 and
            sum(g_v4.GetAtomWithIdx(y).GetSymbol() == "N" for y in r) >= 2 for x in r}
    a, b = np.unravel_index(np.argmax(dm), dm.shape)
    start = min((a, b), key=lambda x: min(dm[x][p] for p in pter))
    order = sorted(by_bond, key=lambda k: min(dm[start][x] for x in k))

    torsions = []
    for n, key in enumerate(order, 1):
        q = by_bond[key]
        torsions.append({"name": f"t{n}", "v4_pos": q, "period": fold_period(g_v4, q),
                         "v4": tuple(hv4[x] for x in q),
                         "mace": tuple(hmc[v4_to_mc[x]] for x in q),
                         "listed_by": "+".join(sorted(source[key]))})
    return torsions, g_v4, v4_sdf


def crystal_values(sdf, g_v4, torsion):
    """Crystal torsion for every symmetry-equivalent mapping onto v4 (folded, unique)."""
    if not os.path.exists(sdf):
        return []
    cry = Chem.MolFromMolFile(sdf, removeHs=False, sanitize=False)
    xyz = cry.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in cry.GetAtoms() if a.GetAtomicNum() > 1]
    g_cry = bare_graph(sdf)
    vals = set()
    for match in g_cry.GetSubstructMatches(g_v4, uniquify=False, maxMatches=1000):
        q = [heavy[match[x]] for x in torsion["v4_pos"]]
        d = dihedral(xyz[None], q)[0]
        vals.add(round(float(fold(d, torsion["period"])), 1))
    return sorted(vals)


def molecule_png(sdf, torsions, width=1400, height=330):
    """2D folate with each torsion's central bond highlighted and labelled t1..tn."""
    mol = Chem.RemoveHs(Chem.MolFromMolFile(sdf))
    AllChem.Compute2DCoords(mol)
    heavy = heavy_index(sdf)
    bonds = []
    for t in torsions:
        j, k = (heavy[x] for x in t["v4_pos"][1:3])
        b = mol.GetBondBetweenAtoms(heavy.index(j), heavy.index(k))
        b.SetProp("bondNote", t["name"])
        bonds.append(b.GetIdx())
    d = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = d.drawOptions()
    opts.annotationFontScale = 0.9
    opts.clearBackground = False
    d.DrawMolecule(mol, highlightAtoms=[], highlightBonds=bonds,
                   highlightBondColors={b: (1.0, 0.85, 0.3) for b in bonds})
    d.FinishDrawing()
    return mpimg.imread(io.BytesIO(d.GetDrawingText()), format="png")


def region(torsion, names):
    """Fragment(s) the central bond belongs to, from v4 atom names."""
    j, k = (names[x] for x in torsion["v4"][1:3])
    hit = [f for f, members in FRAGMENTS.items() if j in members.split() or k in members.split()]
    return "–".join(dict.fromkeys(hit)) or "?"


def plot_system(system, data, torsions, names, crystal, binw, out_png, mol_img):
    ncol = 5
    nrow = int(np.ceil(len(torsions) / ncol))
    fig = plt.figure(figsize=(16, 2.6 + 2.7 * nrow))
    gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[1.0] + [1] * nrow, hspace=0.6, wspace=0.22)
    ax_mol = fig.add_subplot(gs[0, :])
    ax_mol.imshow(mol_img)
    ax_mol.axis("off")
    ns = {lab: data[lab]["ns"] for lab in data}
    fig.suptitle(f"Folate {'in water' if system == 'water' else 'bound to ' + system}, state 0, "
                 f"all frames — v4 (Sage, {ns['v4']:.0f} ns) vs MACE bespoke ({ns['MACE']:.0f} ns)",
                 fontsize=14, y=0.995)
    rows = []
    for n, t in enumerate(torsions):
        ax = fig.add_subplot(gs[1 + n // ncol, n % ncol])
        half = 180.0 / t["period"]
        edges = np.arange(-half, half + binw / 2, binw)
        v4, mc = data["v4"]["tors"][n], data["MACE"]["tors"][n]
        ax.hist(v4, bins=edges, density=True, color=RUNS[0][2], alpha=0.45, lw=0,
                label="v4 (Sage)")
        ax.hist(mc, bins=edges, density=True, histtype="step", color=RUNS[1][2], lw=2,
                label="MACE bespoke")
        for c in crystal[n]:
            ax.axvline(c, color="k", lw=1.2, ls="--")
        ov = overlap(v4, mc, edges)
        self_ov = overlap(v4[: len(v4) // 2], v4[len(v4) // 2:], edges)
        jn, kn = (names[x] for x in t["v4"][1:3])
        fold_note = f"  (mod {360 // t['period']}°)" if t["period"] > 1 else ""
        ax.set_title(f"{t['name']}  {jn}–{kn}  [{region(t, names)}]{fold_note}\n"
                     f"overlap {ov:.2f}   (v4 halves {self_ov:.2f})", fontsize=9.5)
        ax.set_xlim(-half, half)
        ax.set_xticks(np.linspace(-half, half, 5))
        ax.tick_params(labelsize=8)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        rows.append({"system": system, "torsion": t["name"], "central_bond_v4": f"{jn}-{kn}",
                     "region": region(t, names), "period": t["period"],
                     "v4_atoms": " ".join(map(str, t["v4"])), "mace_atoms": " ".join(map(str, t["mace"])),
                     "listed_by": t["listed_by"], "overlap_v4_mace": round(ov, 3),
                     "overlap_v4_halves": round(self_ov, 3),
                     "crystal_deg": " ".join(map(str, crystal[n]))})
    handles, labels = fig.axes[1].get_legend_handles_labels()
    if crystal[0]:
        handles.append(plt.Line2D([], [], color="k", lw=1.2, ls="--"))
        labels.append(f"crystal ({system})")
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.98, 0.97), frameon=False, fontsize=10)
    fig.text(0.5, 0.005, "x: torsion angle (deg); y: density. Dashed black = crystal ligand torsion "
             "(several lines = equivalent atoms). Torsions with equivalent end atoms are folded.",
             ha="center", fontsize=9, color="#555")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--thesis", default=f"{HOME}/thesis-sbrinz-data/thesis")
    ap.add_argument("--crystal", default=f"{HOME}/bsc-thesis/crystal")
    ap.add_argument("--out", default=f"{HOME}/bsc-thesis/results/mace-vs-v4")
    ap.add_argument("--bin", type=float, default=5.0, help="histogram bin width, degrees")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    torsions, g_v4, v4_sdf = build_torsions(args.thesis)
    names = atom_names(f"{args.thesis}/FOL-REST2-v4/water/build/built.solute.pdb")
    mol_img = molecule_png(v4_sdf, torsions)
    for t in torsions:
        print(f"{t['name']:>4}  v4 {t['v4']}  MACE {t['mace']}  period {t['period']}  ({t['listed_by']})")

    rows = []
    for system in SYSTEMS:
        data = {}
        for label, run, _ in RUNS:
            nc = (f"{args.thesis}/{run}/water/REST2-run1/solute_state0_prod1.nc" if system == "water"
                  else f"{args.thesis}/{run}/{system}/analysis-shuyu/FOL_state0_prod1.nc")
            xyz = load_coords(nc)
            f = netcdf_file(nc, "r", mmap=False)
            ns = float(f.variables["time"][-1]) / 1000.0
            f.close()
            key = "v4" if label == "v4" else "mace"
            data[label] = {"ns": ns, "tors": [fold(dihedral(xyz, t[key]), t["period"]) for t in torsions]}
            print(f"{system:>5} {label:>4}: {xyz.shape[0]} frames, {ns:.0f} ns  ({nc})")
        crystal = [crystal_values(f"{args.crystal}/fol_{system}_A.sdf", g_v4, t) if system != "water"
                   else [] for t in torsions]
        png = f"{args.out}/folate_{system}_v4_vs_mace.png"
        rows += plot_system(system, data, torsions, names, crystal, args.bin, png, mol_img)
        print(f"wrote {png}")

    csv_path = f"{args.out}/folate_v4_vs_mace_overlap.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
