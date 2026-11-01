#!/bin/bash
#SBATCH --job-name=rb-inputs
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds/logs/new-inputs-%j.out
# Inputs decided on 2026-11-01, after the first builds:
#   1v2k: the P1 4-piperidinylpyridine as the pyridinium CATION (Open Babel's pH model leaves pyridines
#         neutral). packages/1v2k.input.sdf = the neutral input (kept as 1v2k.input.neutral.sdf) with
#         the pyridine N protonated; Open Babel still runs when the package is made.
#   1l8g (replaces 6cn5) and 6qlp (replaces 5ufr), whose packages md-tools refused: inputs, configs
#         and complex structures with the existing scripts, for these two rows only.
set -euo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
cd "$ROOT"

echo "== 1v2k: pyridinium cation input"
python -I - packages/1v2k.input.neutral.sdf packages/1v2k.input.sdf <<'EOF'
import sys
from rdkit import Chem
m = Chem.MolFromMolFile(sys.argv[1], removeHs=False)
ri = m.GetRingInfo()
pyr = [a for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetIsAromatic() and a.GetDegree() == 2
       and a.GetTotalNumHs() == 0 and ri.IsAtomInRingOfSize(a.GetIdx(), 6)]
assert len(pyr) == 1, f"expected one pyridine N, found {len(pyr)}"
heavy = Chem.RWMol(Chem.RemoveHs(m))
n = heavy.GetAtomWithIdx(pyr[0].GetIdx())        # RemoveHs keeps heavy-atom indices
n.SetFormalCharge(1); n.SetNumExplicitHs(1); n.SetNoImplicit(True)
mol = heavy.GetMol(); Chem.SanitizeMol(mol)
mol = Chem.AddHs(mol, addCoords=True); mol.SetProp("_Name", m.GetProp("_Name"))
Chem.MolToMolFile(mol, sys.argv[2])
print(f"   pyridine N (atom {pyr[0].GetIdx()}) protonated; net charge {Chem.GetFormalCharge(mol)}")
print(f"   {Chem.MolToSmiles(Chem.RemoveHs(mol))}")
EOF

echo "== 1l8g, 6qlp: inputs (CCD chemistry, crystal pose), configs, complex structures"
python -I scripts/make_inputs.py "$ROOT" 1l8g 6qlp
python -I scripts/make_configs.py "$ROOT" 1l8g 6qlp
python -I scripts/prep.py "$ROOT" 1l8g 6qlp
echo "== done $(date)"
