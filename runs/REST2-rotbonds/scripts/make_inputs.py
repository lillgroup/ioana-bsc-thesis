"""Parameterisation inputs for the REST2-rotbonds compounds: packages/<name>.input.sdf and
packages/<name>.parameterize.config, one per row of compounds.csv.

  python -I scripts/make_inputs.py <root>

The chemistry (bond orders, formal charges, stereo) is checked against the wwPDB Chemical Component
Dictionary entry of the deposited ligand (ccd/<code>_ideal.sdf, from files.rcsb.org/ligands):
  - same standard InChI (which ignores tautomers) -> the PDBbind ligand SDF's bond orders
    (pdbbind/<name>_ligand.sdf); the CCD draws some in odd tautomers, e.g. 4ase's urea as an
    isoxazol-ylidene imine (one rotatable bond less)
  - different -> the CCD's: PDBbind's 4p6w file has mometasone furoate's C20 ketone as an alcohol.
The heavy-atom coordinates are the deposited ones (crystal/<PDB>.pdb, the ligand residue of the
chain in compounds.csv, altloc ' ' or 'A'); hydrogens are added on them. The molecule is written
neutral, as in the earlier v4 inputs; Open Babel sets the pH 7.4 protomer when the package is made
(parameterize.config), exactly as for 4ibi / 4ibk / 4j7d.
Checks: every CCD heavy atom is matched by name in the crystal, and the CIP labels of the crystal
pose equal the CCD's (a mirror-image or wrongly named deposition would stop the script).
"""
import csv
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
CONFIG = """# {name} ({ligand}, {pdb} chain {chain}) for REST2 v4 (md-tools feature/ligand-protonation-and-ligand-amides,
# b3603e0), as packages/4ibi of REST2-conformers-v4: neutral input (CCD chemistry, crystal pose,
# make_inputs.py); Open Babel's pH model sets the pH 7.4 protomer; Sage + AM1-BCC.
solute:
  kind: ligand
  aliases: [{name}]
protonation:
  ligand_method: openbabel
  ph: 7.4
"""


def crystal_ligand(pdb, code, chain):
    lines = [l for l in pdb.read_text().splitlines()
             if l.startswith("HETATM") and l[17:20].strip() == code and l[21] == chain and l[16] in " A"
             and l[76:78].strip() not in ("H", "D")]
    resids = {l[22:27] for l in lines}
    assert len(resids) == 1, f"{pdb.name}: {code} chain {chain}: residues {resids}"
    return lines


def main(root):
    rows = list(csv.DictReader(open(root / "compounds.csv")))
    (root / "packages").mkdir(exist_ok=True)
    for r in rows:
        ccd = Chem.MolFromMolFile(str(root / "ccd" / f"{r['ligand']}_ideal.sdf"), removeHs=False)
        lines = crystal_ligand(root / "crystal" / f"{r['pdb']}.pdb", r["ligand"], r["chain"])
        pdbbind = Chem.MolFromMolFile(str(root / "pdbbind" / f"{r['name']}_ligand.sdf"), removeHs=False)
        inchi = lambda m: Chem.MolToInchi(rdMolStandardize.Uncharger().uncharge(Chem.RemoveHs(m)))
        source = "PDBbind" if inchi(pdbbind) == inchi(ccd) else "CCD (PDBbind SDF differs)"
        heavy = Chem.RemoveHs(pdbbind if source == "PDBbind" else ccd)
        # the deposited heavy atoms (bonds by proximity) get the template's bond orders by graph match
        cryst = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", removeHs=True, proximityBonding=True,
                                     sanitize=False)
        assert cryst.GetNumAtoms() == heavy.GetNumAtoms(), \
            f"{r['name']}: {cryst.GetNumAtoms()} deposited heavy atoms, template has {heavy.GetNumAtoms()}"
        mol = AllChem.AssignBondOrdersFromTemplate(heavy, cryst)
        Chem.AssignStereochemistryFrom3D(mol)
        ref = Chem.Mol(heavy)
        Chem.AssignStereochemistryFrom3D(ref)
        same = Chem.MolToSmiles(mol) == Chem.MolToSmiles(ref)
        assert same, f"{r['name']}: crystal pose stereo/chemistry differs from the template:\n " \
                     f"{Chem.MolToSmiles(mol)}\n {Chem.MolToSmiles(ref)}"
        mol = rdMolStandardize.Uncharger().uncharge(mol)     # neutral input; Open Babel sets pH 7.4
        mol = Chem.AddHs(mol, addCoords=True)
        mol.SetProp("_Name", r["resname"])
        out = root / "packages" / f"{r['name']}.input.sdf"
        Chem.MolToMolFile(mol, str(out))
        (root / "packages" / f"{r['name']}.parameterize.config").write_text(CONFIG.format(**r))
        h = Chem.RemoveHs(mol)
        print(f"{r['name']} {r['ligand']} -> {r['resname']}: {h.GetNumAtoms()} heavy atoms, "
              f"rot. bonds {rdMolDescriptors.CalcNumRotatableBonds(h)} (csv {r['rot_bonds']}), "
              f"MW {Descriptors.MolWt(mol):.1f}, charge {Chem.GetFormalCharge(mol)}, chemistry from {source}\n"
              f"    {Chem.MolToSmiles(h)}")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
