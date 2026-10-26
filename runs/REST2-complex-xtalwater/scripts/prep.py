"""Complex inputs for REST2, from the deposited crystal structures, crystal waters included.

  python -I scripts/prep.py <root> [name ...]

For each system, <root>/<name>/prep/complex.pdb holds, from <root>/crystal/<PDB>.pdb:
  chain A protein   (altloc ' ' or 'A', altloc column blanked; biological assembly 1 is the chain-A
                     monomer for every entry below)
  TER
  chain-A ligand(s) heavy atoms as deposited (md-tools adds the package hydrogens)
  TER
  chain-A crystal waters (oxygens; md-tools adds hydrogens, keeps them after the solute and fills
                     the rest of the box around them)
Every other HETATM group of chain A must be listed in EXCLUDED with its reason, so nothing is
dropped silently; an unlisted one stops the script. The exclusions were approved by Shu-Yu on
2026-10-26 (README.md). Other chains are not read at all.

This replaces bsc-thesis/runs/FOL-complex-REST2/scripts/prep.py and
REST2-conformers-complex/scripts/prep.py, which dropped the crystal waters. The protein and ligand
records written here are identical to theirs.
"""
import math
import sys
from pathlib import Path

# name: (PDB id, ligand residue names kept)
SYSTEMS = {
    "FOL-4LRH": ("4LRH", ["FOL"]),
    "FOL-3QL3": ("3QL3", ["FOL", "NAP"]),
    "4ibi": ("4IBI", ["1D8"]),
    "4ibk": ("4IBK", ["1DE"]),
    "4j7d": ("4J7D", ["I31"]),
    "5mwh": ("5MWH", ["UWX"]),
    "5ov8": ("5OV8", ["AXN"]),
}
# chain-A HETATM groups left out, and why
EXCLUDED = {
    "DMS": "crystallisation solution (DMSO)",
    "GOL": "crystallisation solution (glycerol)",
    "SO4": "crystallisation solution (sulfate)",
    "NO3": "crystallisation solution (nitrate)",
    "MN": "force field: no Mn2+ parameters in md-tools; far from the ligand site",
    "NAG": "force field: glycosylated Asn built as Asn (no carbohydrate parameters); far from the site",
}


def xyz(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def prepare(root, name):
    pdb, ligands = SYSTEMS[name]
    prot, lig, water, left_out = [], [], [], {}
    for line in (root / "crystal" / f"{pdb}.pdb").read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")) or line[21] != "A" or line[16] not in " A":
            continue
        line = line[:16] + " " + line[17:]
        res = line[17:20].strip()
        if line.startswith("ATOM"):
            prot.append(line)
        elif res in ligands:
            lig.append(line)
        elif res == "HOH":
            water.append(line)
        elif res in EXCLUDED:
            left_out.setdefault((res, line[22:27].strip()), []).append(xyz(line))
        else:
            sys.exit(f"{name}: {pdb} chain A has {res} {line[22:27].strip()}, which is neither a "
                     f"ligand nor listed in EXCLUDED; decide what to do with it first")
    assert {l[17:20].strip() for l in lig} == set(ligands), f"{name}: ligands missing"
    assert all(l[12:16].strip() == "O" for l in water), f"{name}: expected oxygen-only waters"

    # A real overlap, not a short crystal hydrogen bond (several deposited waters sit 2.15-2.4 A from
    # Asp/Arg/Asn N/O atoms): closer than 1.8 A to anything, or 2.6 A to a carbon.
    solute = [(xyz(l), l[76:78].strip()) for l in prot + lig]
    clash = [w for w in water
             if any(math.dist(xyz(w), p) < (2.6 if el == "C" else 1.8) for p, el in solute)]
    assert not clash, f"{name}: waters overlapping the solute: {[w[17:27] for w in clash]}"

    lig_xyz = [xyz(l) for l in lig]
    near = sum(min(math.dist(xyz(w), p) for p in lig_xyz) < 3.5 for w in water)
    print(f"{name} ({pdb}): {len(prot)} protein atoms, {len(lig)} ligand atoms "
          f"({' '.join(ligands)}), {len(water)} crystal waters ({near} within 3.5 A of the ligands)")
    for (res, resid), pts in sorted(left_out.items()):
        d = min(math.dist(p, q) for p in pts for q in lig_xyz)
        print(f"  left out: {res} A{resid}, {d:.1f} A from the ligands -- {EXCLUDED[res]}")

    out = root / name / "prep" / "complex.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    remark = (f"REMARK   1 {pdb} chain A + {' + '.join(ligands)} + {len(water)} crystal waters; "
              f"left out: {', '.join(sorted({r for r, _ in left_out})) or 'nothing'} (prep.py)")
    out.write_text("\n".join([remark] + prot + ["TER"] + lig + ["TER"] + water + ["END"]) + "\n")


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve()
    for name in sys.argv[2:] or SYSTEMS:
        prepare(root, name)
