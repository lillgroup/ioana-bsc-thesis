"""Add the deposited crystal waters to a complex prepared without them.

  python -I scripts/add_crystal_waters.py <system dir> <crystal .pdb> [--chain A]

The earlier complex runs started from prep/complex.pdb written by prep.py, which kept chain A and
its ligand(s) but dropped every crystal water, so the solvation step refilled the pocket with bulk
water. md-tools' own protocol keeps the deposited waters (they are moved after the solute and
`addSolvent` fills the rest; see md_tools.openmm.solvation._retained_solvent_last). This writes

  prep/complex.pdb = prep/complex.nowater.pdb (protein + ligands, unchanged, line for line)
                     + every HOH of the chain (altloc ' ' or 'A', altloc column blanked)

Nothing else is added back: ions md-tools has no parameters for (3QL3 Mn2+) and crystallisation
additives (DMS, GOL, SO4, NO3) stay out, as do glycans (4LRH NAG); their distances to the ligand
are printed for the record. Every protein/ligand atom of complex.nowater.pdb is first checked
against the deposit (same residue, atom name, coordinates), so the waters are placed in the frame
the complex was taken from.
"""
import argparse
import math
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("system")
ap.add_argument("crystal")
ap.add_argument("--chain", default="A")
args = ap.parse_args()
prep = Path(args.system) / "prep"
src = (prep / "complex.nowater.pdb").read_text().splitlines()


def key(line):
    return line[17:20].strip(), line[21], line[22:27].strip(), line[12:16].strip()


def xyz(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


crystal = [l for l in Path(args.crystal).read_text().splitlines()
           if l.startswith(("ATOM", "HETATM")) and l[16] in " A"]
by_key = {key(l): xyz(l) for l in crystal}

atoms = [l for l in src if l.startswith(("ATOM", "HETATM"))]
missing = [l for l in atoms if key(l) not in by_key]
moved = [l for l in atoms if key(l) in by_key and math.dist(xyz(l), by_key[key(l)]) > 1e-3]
assert not missing and not moved, (f"complex.nowater.pdb does not match {args.crystal}: "
                                   f"{len(missing)} atoms missing, {len(moved)} moved")
print(f"{len(atoms)} protein/ligand atoms match {Path(args.crystal).name} exactly")

ligands = {l[17:20] for l in atoms if l.startswith("HETATM")}
lig_xyz = [xyz(l) for l in atoms if l.startswith("HETATM")]
prot_xyz = [xyz(l) for l in atoms if l.startswith("ATOM")]
waters = [l[:16] + " " + l[17:] for l in crystal
          if l.startswith("HETATM") and l[17:20] == "HOH" and l[21] == args.chain]
assert waters and all(l[12:16].strip() == "O" for l in waters), "expected oxygen-only waters"
near = sum(min(math.dist(xyz(w), p) for p in lig_xyz) < 3.5 for w in waters)
# A real overlap, not a short crystal hydrogen bond (5OV8 has three deposited waters 2.15-2.18 A
# from Asp/Arg/Asn N/O atoms): closer than 1.8 A to anything, or 2.6 A to a carbon.
solute = [(xyz(l), l[76:78].strip()) for l in atoms]
clash = [w for w in waters
         if any(math.dist(xyz(w), p) < (2.6 if el == "C" else 1.8) for p, el in solute)]
assert not clash, f"waters overlapping the solute: {[key(w) for w in clash]}"
short = [w for w in waters if min(math.dist(xyz(w), p) for p, _ in solute) < 2.4]
for w in short:
    d, el = min((math.dist(xyz(w), p), el) for p, el in solute)
    print(f"  kept: {' '.join(key(w)[:3])} {d:.2f} A from a solute {el} (short hydrogen bond)")
print(f"adding {len(waters)} chain-{args.chain} crystal waters ({near} within 3.5 A of the ligand(s))")

dropped = {}
for l in crystal:
    if l.startswith("HETATM") and l[17:20] != "HOH" and l[17:20] not in ligands:
        dropped.setdefault((l[17:20].strip(), l[21], l[22:27].strip()), []).append(xyz(l))
for (res, chain, resid), pts in sorted(dropped.items()):
    d = min(math.dist(p, q) for p in pts for q in lig_xyz)
    print(f"  not added: {res} {chain}{resid}  {d:.1f} A from the ligand(s)")

body = [l for l in src if not l.startswith("END")]
remark = (f"REMARK   1 + {len(waters)} chain-{args.chain} crystal waters from "
          f"{Path(args.crystal).name} (add_crystal_waters.py)")
out = [body[0], remark] + body[1:] + ["TER"] + waters + ["END"]
(prep / "complex.pdb").write_text("\n".join(out) + "\n")
print(f"wrote {prep / 'complex.pdb'}")
