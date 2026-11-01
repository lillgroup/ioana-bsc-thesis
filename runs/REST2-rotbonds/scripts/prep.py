"""Complex inputs for REST2-rotbonds, from the deposited crystal structures, crystal waters included.
Same protocol as REST2-complex-xtalwater/scripts/prep.py, for the chain named in compounds.csv.

  python -I scripts/prep.py <root> [name ...]

For each compound, <root>/<name>-complex/prep/complex.pdb holds, from <root>/crystal/<PDB>.pdb:
  the ligand's chain, protein (altloc ' ' or 'A', altloc column blanked; deposited H/D dropped)
  TER
  the ligand heavy atoms as deposited (md-tools adds the package hydrogens)
  TER
  that chain's crystal waters (oxygens) and ions (md-tools keeps both, after the solute, and fills
                     the rest of the box around them)
Every other HETATM group of the chain must be in EXCLUDED with its reason, so nothing is dropped
silently; an unlisted one stops the script. Other chains are not read (biological assembly 1 is
one chain for all 15; checked in the selection: no other chain within 5 A of a ligand is a bound
peptide, only crystal packing). md-tools refuses residues missing INSIDE a chain, so the selection
took only structures without such gaps; this script checks again.
"""
import csv
import math
import sys
from pathlib import Path

ION_NAMES = {"NA", "CL", "K", "MG", "CA", "ZN", "BR", "I"}   # md-tools keeps deposited ions
# HETATM groups left out, and why (crystallisation / cryo-protection additives, as in
# REST2-complex-xtalwater, where Shu-Yu approved leaving them out on 2026-10-26)
EXCLUDED = {
    "SO4": "crystallisation solution (sulfate)",
    "PO4": "crystallisation solution (phosphate)",
    "GOL": "cryo-protectant (glycerol)",
    "EDO": "cryo-protectant (ethylene glycol)",
    "PG4": "crystallisation solution (PEG fragment)",
    "ACY": "crystallisation solution (acetate)",
    "DMS": "crystallisation solution (DMSO)",
}

# single deposited waters left out, and why (the clash check below stops on any other)
DROP_WATERS = {
    ("1eld", "651"): "2.06 A from ligand C5 (and 2.6 A from C4/C6): overlaps the inhibitor, a deposition error",
}


def xyz(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def chain_breaks(prot):
    res = {}
    for l in prot:
        res.setdefault((l[22:26], l[26]), {})[l[12:16].strip()] = xyz(l)
    keys = list(res)
    return [f"{a[0].strip()}{a[1].strip()}-{b[0].strip()}{b[1].strip()}" for a, b in zip(keys, keys[1:])
            if "C" in res[a] and "N" in res[b] and math.dist(res[a]["C"], res[b]["N"]) > 2.0]


def prepare(root, c):
    name, pdb, ligand, chain = c["name"], c["pdb"], c["ligand"], c["chain"]
    prot, lig, water, ions, left_out = [], [], [], [], {}
    for line in (root / "crystal" / f"{pdb}.pdb").read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")) or line[21] != chain or line[16] not in " A":
            continue
        if line[76:78].strip() in ("H", "D"):
            continue
        line = line[:16] + " " + line[17:]
        res = line[17:20].strip()
        if line.startswith("ATOM"):
            prot.append(line)
        elif res == ligand:
            lig.append(line)
        elif res == "HOH" and (name, line[22:27].strip()) in DROP_WATERS:
            print(f"  left out: HOH {chain}{line[22:27].strip()} -- {DROP_WATERS[(name, line[22:27].strip())]}")
        elif res == "HOH":
            water.append(line)
        elif res in ION_NAMES:
            ions.append(line)
        elif res in EXCLUDED:
            left_out.setdefault((res, line[22:27].strip()), []).append(xyz(line))
        else:
            sys.exit(f"{name}: {pdb} chain {chain} has {res} {line[22:27].strip()}, which is neither the "
                     f"ligand nor listed in EXCLUDED; decide what to do with it first")
    assert lig, f"{name}: no {ligand} in chain {chain}"
    assert len({l[22:27] for l in lig}) == 1, f"{name}: several {ligand} residues in chain {chain}"
    assert all(l[12:16].strip() == "O" for l in water), f"{name}: expected oxygen-only waters"
    breaks = chain_breaks(prot)
    assert not breaks, f"{name}: residues missing inside chain {chain} ({breaks}); md-tools refuses gaps"

    solute = [(xyz(l), l[76:78].strip()) for l in prot + lig]
    clash = [w for w in water
             if any(math.dist(xyz(w), p) < (2.6 if el == "C" else 1.8) for p, el in solute)]
    assert not clash, f"{name}: waters overlapping the solute: {[w[17:27] for w in clash]}"

    lig_xyz = [xyz(l) for l in lig]
    near = sum(min(math.dist(xyz(w), p) for p in lig_xyz) < 3.5 for w in water)
    print(f"{name} ({pdb} chain {chain}): {len({l[22:27] for l in prot})} residues, {len(lig)} ligand atoms "
          f"({ligand}), {len(water)} crystal waters ({near} within 3.5 A of the ligand), "
          f"ions: {' '.join(sorted(l[17:20].strip() for l in ions)) or 'none'}")
    for (res, resid), pts in sorted(left_out.items()):
        d = min(math.dist(p, q) for p in pts for q in lig_xyz)
        print(f"  left out: {res} {chain}{resid}, {d:.1f} A from the ligand -- {EXCLUDED[res]}")

    out = root / f"{name}-complex" / "prep" / "complex.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    remark = (f"REMARK   1 {pdb} chain {chain} + {ligand} + {len(water)} crystal waters + {len(ions)} ions; "
              f"left out: {', '.join(sorted({r for r, _ in left_out})) or 'nothing'} (prep.py)")
    out.write_text("\n".join([remark] + prot + ["TER"] + lig + ["TER"] + water + ions + ["END"]) + "\n")


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve()
    rows = {r["name"]: r for r in csv.DictReader(open(root / "compounds.csv"))}
    for name in sys.argv[2:] or rows:
        prepare(root, rows[name])
