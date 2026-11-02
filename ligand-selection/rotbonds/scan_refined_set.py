"""Step 1 of the REST2-rotbonds selection: properties of every PDBbind v2020 refined-set complex.

  python -I ligand-selection/rotbonds/scan_refined_set.py <bsc-thesis>/refined-set ligand-selection/rotbonds/refined_set_properties.csv

Per complex: ligand MW, rotatable bonds (rb = filter.py's definition, on the SDF WITH explicit H, which
also counts OH/NH rotors; rb_heavy = the usual one, on the heavy-atom graph), heavy atoms, formal
charge, rings, elements; protein residues/chains/atoms; non-water HETATM groups in the pocket file
(metals flagged); UniProt id, name, resolution and affinity from the index files. The refined set
itself is not in git (.gitignore; PDBbind licence) -- it is on the local machine only.
"""
import os, sys
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')
R = sys.argv[1]
name = {}
for l in open(f"{R}/index/INDEX_refined_name.2020"):
    if l.startswith("#"): continue
    p = l.split(None, 3); name[p[0]] = (p[2], p[3].strip())
data = {}
for l in open(f"{R}/index/INDEX_refined_data.2020"):
    if l.startswith("#"): continue
    p = l.split(); data[p[0]] = (p[1], p[4], l.split("(")[-1].rstrip(")\n") if "(" in l else "")
METALS = {"ZN","MG","CA","MN","FE","FE2","CO","CU","NI","NA","K","CD","HG"}
rows = []
for code in sorted(os.listdir(R)):
    sdf = f"{R}/{code}/{code}_ligand.sdf"
    if not os.path.isfile(sdf): continue
    sup = Chem.SDMolSupplier(sdf, removeHs=False)
    mol = sup[0] if len(sup) else None
    if mol is None: continue
    heavy = Chem.RemoveHs(mol)
    res, chains, atoms = set(), set(), 0
    for l in open(f"{R}/{code}/{code}_protein.pdb"):
        if l.startswith("ATOM"):
            atoms += 1
            if l[12:16].strip() == "CA": res.add((l[21], l[22:27])); chains.add(l[21])
    het = set()
    for l in open(f"{R}/{code}/{code}_pocket.pdb"):
        if l.startswith("HETATM"):
            rn = l[17:20].strip()
            if rn != "HOH": het.add(rn)
    rows.append(dict(pdb=code, mw=Descriptors.MolWt(mol),
        rb=rdMolDescriptors.CalcNumRotatableBonds(mol),            # filter.py's definition (explicit H)
        rb_heavy=rdMolDescriptors.CalcNumRotatableBonds(heavy),
        ha=mol.GetNumHeavyAtoms(), charge=Chem.GetFormalCharge(mol),
        rings=rdMolDescriptors.CalcNumRings(mol),
        elements="".join(sorted({a.GetSymbol() for a in heavy.GetAtoms()})),
        n_res=len(res), n_chains=len(chains), n_atoms=atoms,
        pocket_het=" ".join(sorted(het)), pocket_metal=bool(het & METALS),
        uniprot=name.get(code, ("", ""))[0], protein=name.get(code, ("", ""))[1],
        resolution=data.get(code, ("",))[0], affinity=data.get(code, ("", "", ""))[1],
        ligname=data.get(code, ("", "", ""))[2], smiles=Chem.MolToSmiles(heavy)))
pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
print(len(rows))
