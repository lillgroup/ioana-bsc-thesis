"""Inputs for one compound's report, in <work>/<compound>/ (outside git):

  water.nc / water.pdb / water.sdf     state-0 ligand trajectory of the water run (run + extensions
                                       joined), its topology and the package SDF (bond orders)
  <label>.nc / <label>.pdb / .sdf      the same for each complex: the ligand cut out of the state-0
                                       solute trajectory, in the complex package's own atom order
  torsions/<system>/torsions{,_rigid,_polarH}.txt
                                       ONE torsion definition for every plot of the compound: written
                                       for the water ligand by select_torsions.py (rotatable, rigid,
                                       polar-H; ring torsions are not used), then mapped onto each
                                       complex ligand by graph matching -- so t<k> is the same bond in
                                       every system, also when the complex package orders its atoms
                                       differently (3QL3: folate 3H tautomer, a separate package)
  manifest.json                        everything the report steps read

  python -I prepare.py --name <compound> [--table compounds.csv] [--thesis DIR] [--repo DIR] [--work DIR]

Paths in compounds.csv: run folders relative to --thesis; crystal references relative to --repo.
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys

import mdtraj as md
from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
SHUYU = os.path.join(os.path.dirname(HERE), "shuyu_clustering")
KINDS = ("torsions", "torsions_rigid", "torsions_polarH")


def ligand_files(run):
    """(residue name, package SDF) of a run: the build's one SDF whose stem is a residue of
    built.solute.pdb (a water build may also hold the parameterisation input, e.g. 5mwh.sdf)."""
    residues = {r.name for r in md.load(f"{run}/build/built.solute.pdb").topology.residues}
    hits = [p for p in glob.glob(f"{run}/build/*.sdf") if os.path.basename(p)[:-4] in residues]
    if len(hits) != 1:
        sys.exit(f"{run}/build: expected one SDF named after a solute residue, found {hits}")
    return os.path.basename(hits[0])[:-4], hits[0]


def skeleton(mol):
    """Heavy-atom graph, all bonds single, no charges or aromaticity: identical for tautomers and
    protomers of one molecule. Returns (graph, heavy atom indices of mol)."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rw, new = Chem.RWMol(), {}
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


def atom_map(water, cplx):
    """water atom index -> complex atom index (None where no counterpart, e.g. a tautomeric H).
    Heavy atoms by graph match (the identity when both are the same package); H by parent atom and
    rank among the parent's H."""
    gw, hw = skeleton(water)
    gc, hc = skeleton(cplx)
    if gw.GetNumAtoms() != gc.GetNumAtoms():
        sys.exit(f"heavy atoms differ: water {gw.GetNumAtoms()}, complex {gc.GetNumAtoms()}")
    matches = gc.GetSubstructMatches(gw, uniquify=False, useChirality=False, maxMatches=100000)
    if not matches:
        sys.exit("the complex ligand is not the water ligand (heavy-atom graphs differ)")
    ident = [m for m in matches if all(hc[m[i]] == hw[i] for i in range(len(hw)))]
    m = ident[0] if ident else matches[0]
    amap = {hw[i]: hc[m[i]] for i in range(len(hw))}
    hs = lambda mol, i: sorted(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1)
    for i in hw:
        for r, h in enumerate(hs(water, i)):
            ch = hs(cplx, amap[i])
            amap[h] = ch[r] if r < len(ch) else None
    return amap, bool(ident)


def read_quads(path):
    out = []
    for line in open(path):
        if line.strip() and not line.strip().startswith("#"):
            body, _, comment = line.partition("#")
            out.append((tuple(int(x) for x in body.split()[:4]), comment.strip()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--name", required=True)
    ap.add_argument("--table", default=os.path.join(HERE, "compounds.csv"))
    ap.add_argument("--thesis", default="/scicore/home/lill0000/popa0002/thesis")
    ap.add_argument("--repo", default=os.path.dirname(os.path.dirname(HERE)))
    ap.add_argument("--work", default="/scicore/home/lill0000/popa0002/thesis/analysis-report")
    args = ap.parse_args()
    row = {r["name"]: r for r in csv.DictReader(open(args.table))}[args.name]
    work = os.path.join(args.work, args.name)
    os.makedirs(work, exist_ok=True)
    T = lambda p: os.path.join(args.thesis, p)
    R = lambda p: os.path.join(args.repo, p)

    # ---- water: trajectory (run + extensions), topology, SDF, torsions -------------------------
    wrun = T(row["water"])
    wres, wsdf = ligand_files(wrun)
    parts = [f"{wrun}/REST2-run1/solute_state0_prod1.nc"] + \
            [f"{wrun}/{e}/solute_state0_prod1.nc" for e in filter(None, row["water_ext"].split(";"))]
    t = md.join([md.load(p, top=f"{wrun}/build/built.solute.pdb") for p in parts])
    t.save_netcdf(f"{work}/water.nc")
    t[0].save_pdb(f"{work}/water.pdb")
    subprocess.run(["cp", wsdf, f"{work}/water.sdf"], check=True)
    print(f"water: {row['water']} ({', '.join(os.path.relpath(p, wrun) for p in parts)}), "
          f"{t.n_frames} frames x {t.n_atoms} atoms, ligand {wres}")
    tdir = f"{work}/torsions/water"
    os.makedirs(tdir, exist_ok=True)
    subprocess.run([sys.executable, "-I", f"{SHUYU}/select_torsions.py", "--sdf", wsdf,
                    "--top", f"{wrun}/build/built.pdb", "--select", f"resname {wres}", "--ligand-indices",
                    "--polar-h", "--out", f"{tdir}/torsions"], check=True, stdout=subprocess.DEVNULL)
    wquads = {k: read_quads(f"{tdir}/{k}.txt") for k in KINDS}
    numbering = [q for k in KINDS for q, _ in wquads[k]]
    print(f"torsions: {', '.join(f'{len(wquads[k])} {k}' for k in KINDS)} -> t0..t{len(numbering) - 1}")
    water_mol = Chem.MolFromMolFile(wsdf, removeHs=False, sanitize=False)
    assert [a.element.symbol.upper() for a in t.topology.atoms] == \
           [a.GetSymbol().upper() for a in water_mol.GetAtoms()], "water topology and SDF atom order differ"

    manifest = {"name": args.name, "version": row["version"], "stride": int(row["stride"]),
                "max_matches": row["max_matches"] or None,
                "water": {"label": "water", "run": wrun, "traj": f"{work}/water.nc", "top": f"{work}/water.pdb",
                          "sdf": f"{work}/water.sdf", "torsions": tdir,
                          "refs": [R(p) for p in filter(None, row["water_refs"].split(";"))],
                          "rest2_config": f"{wrun}/REST2.config" if os.path.exists(f"{wrun}/REST2.config") else None,
                          "parts": parts},
                "complexes": []}

    # ---- complexes: cut the ligand, map the torsions ---------------------------------------------
    for spec in filter(None, row["complexes"].split(";")):
        label, run, ref = spec.split("=")
        crun = T(run)
        cres, csdf = ligand_files(crun)
        full = md.load(f"{crun}/build/built.solute.pdb")
        # by name in Python: mdtraj's selection language cannot parse residue names starting with a digit (1D8)
        lig = [a.index for a in full.topology.atoms if a.residue.name == cres]
        cmol = Chem.MolFromMolFile(csdf, removeHs=False, sanitize=False)
        if len(lig) != cmol.GetNumAtoms():
            sys.exit(f"{label}: {len(lig)} atoms of residue {cres}, package SDF has {cmol.GetNumAtoms()}")
        assert [full.topology.atom(i).element.symbol.upper() for i in lig] == \
               [a.GetSymbol().upper() for a in cmol.GetAtoms()], f"{label}: ligand atom order differs from its SDF"
        c = md.load(f"{crun}/REST2-run1/solute_state0_prod1.nc", top=f"{crun}/build/built.solute.pdb",
                    atom_indices=lig)
        c.save_netcdf(f"{work}/{label}.nc")
        for a in c.topology.atoms:     # renumber: the cut keeps the complex's serials (e.g. 1938..), the
            a.serial = None            # CONECT records count from 1, and the bonds would be lost on reading
        c[0].save_pdb(f"{work}/{label}.pdb")
        subprocess.run(["cp", csdf, f"{work}/{label}.sdf"], check=True)
        amap, identity = atom_map(Chem.MolFromMolFile(wsdf, removeHs=False), Chem.MolFromMolFile(csdf, removeHs=False))
        cdir = f"{work}/torsions/{label}"
        os.makedirs(cdir, exist_ok=True)
        for k in KINDS:
            lines = []
            for q, comment in wquads[k]:
                mq = [amap.get(i) for i in q]
                if None in mq:
                    sys.exit(f"{label}: torsion {q} ({k}) uses an atom with no counterpart in the complex ligand")
                lines.append(" ".join(map(str, mq)) + f"  # water {' '.join(map(str, q))}" + (f"  {comment}" if comment else ""))
            open(f"{cdir}/{k}.txt", "w").write("\n".join(lines) + ("\n" if lines else ""))
        print(f"{label}: {run}, {c.n_frames} frames x {c.n_atoms} atoms ({cres}); torsions mapped "
              f"{'1:1 (same package and atom order)' if identity else 'by graph matching (different atom order)'}")
        manifest["complexes"].append({"label": label, "run": crun, "run_rel": run, "traj": f"{work}/{label}.nc",
                                      "top": f"{work}/{label}.pdb", "sdf": f"{work}/{label}.sdf",
                                      "torsions": cdir, "refs": [R(ref)], "identity_map": identity,
                                      "rest2_config": f"{crun}/REST2.config"})
    json.dump(manifest, open(f"{work}/manifest.json", "w"), indent=2)
    print(f"wrote {work}/manifest.json")


if __name__ == "__main__":
    main()
