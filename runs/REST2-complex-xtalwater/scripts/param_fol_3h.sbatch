#!/bin/bash
#SBATCH --job-name=xw-param-3H
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/param-fol-3h-%j.out
# Folate parameter package for the 3QL3 TAUTOMER TEST: the v4 folate with the pteridine H on N3
# instead of N1 (crystal: Asp27 OD2 2.64 A from N3, so N3 must donate). Same route as the v4 package:
# neutral SMILES -> Open Babel pH 7.4 (both COO-, net -2) -> Sage 2.2.1 + AM1-BCC. New species, so a
# new package (approved by the user 2026-10-29); the v4 package FOL is not touched.
# Checks: (1) the input differs from the v4 input only in where the H sits (same heavy-atom graph and
# glutamate stereo); (2) after Open Babel the package is still H on N3, none on N1, net -2.
set -uo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
V4IN=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/runs/FOL-REST2-v4/packages/FOL.input.sdf
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
echo "== $(md-openmm --version)  on $(hostname)  $(date)"
fail() { echo "FAILED: $1"; exit 1; }
cd "$ROOT/packages"
[ -e FOL-3H ] && fail "packages/FOL-3H exists already; not overwriting it"

echo "== 1. input = v4 folate input with the H moved from N1 to N3"
python -I - FOL-3H.smi "$V4IN" <<'EOF' || fail "input check"
import sys
from rdkit import Chem

def skeleton(m):
    """Heavy-atom graph (all bonds single, neutral, no H) + CIP labels: identical for tautomers."""
    m = Chem.RemoveHs(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    cip = sorted(a.GetProp("_CIPCode") for a in m.GetAtoms() if a.HasProp("_CIPCode"))
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False); a.SetFormalCharge(0); a.SetNoImplicit(True); a.SetNumExplicitHs(0)
    return Chem.MolToSmiles(rw), cip

new = Chem.MolFromSmiles(open(sys.argv[1]).read().split()[0])
old = Chem.MolFromMolFile(sys.argv[2])
same = skeleton(new) == skeleton(old)
print(f"   heavy-atom graph and stereo identical to the v4 input: {same}  (CIP {skeleton(new)[1]})")
print(f"   standard InChI identical (it ignores this tautomerism): {Chem.MolToInchi(new) == Chem.MolToInchi(old)}")
sys.exit(0 if same else 1)
EOF

echo "== 2. parameterise (Open Babel pH 7.4, AM1-BCC via sqm)"
md-openmm build-top --parameterize -i FOL-3H.smi --resname FOL --config FOL-3H.parameterize.config \
    -op FOL-3H/FOL.pdb -os FOL-3H/FOL.xml -log FOL-3H.parameterize.log > FOL-3H.parameterize.out 2>&1 \
  || { grep -v '^$' FOL-3H.parameterize.out | tail -20; fail parameterise; }
grep -E 'net_formal_charge|canonical_smiles|parameter_id' FOL-3H/parameter.config
grep -iE 'protonation' FOL-3H.parameterize.out | head -3

echo "== 3. the package is the 3H dianion (and v4 is the N1-H dianion)"
python -I - FOL-3H/molecule.sdf FOL/molecule.sdf <<'EOF' || fail "the new package is not the 3H dianion"
import sys
from rdkit import Chem

def tautomer(path):
    """(H on N1, H on N3, net charge). C2 = ring C with an exocyclic N and two ring N; N3 is the
    ring N next to the C=O carbon (C4), N1 the other one."""
    m = Chem.MolFromMolFile(path, removeHs=False)
    nh = lambda a: sum(x.GetSymbol() == "H" for x in a.GetNeighbors())
    def next_to_carbonyl(n, c2):
        for c in n.GetNeighbors():
            if c.GetIdx() != c2.GetIdx() and c.GetSymbol() == "C":
                for o in c.GetNeighbors():
                    b = m.GetBondBetweenAtoms(c.GetIdx(), o.GetIdx())
                    if o.GetSymbol() == "O" and b.GetBondTypeAsDouble() == 2:
                        return True
        return False
    for c2 in m.GetAtoms():
        ring_n = [x for x in c2.GetNeighbors() if x.GetSymbol() == "N" and x.IsInRing()]
        exo_n = [x for x in c2.GetNeighbors() if x.GetSymbol() == "N" and not x.IsInRing()]
        if c2.GetSymbol() == "C" and c2.IsInRing() and len(ring_n) == 2 and len(exo_n) == 1:
            n3 = [x for x in ring_n if next_to_carbonyl(x, c2)]
            if len(n3) == 1:
                n1 = [x for x in ring_n if x.GetIdx() != n3[0].GetIdx()][0]
                return nh(n1), nh(n3[0]), Chem.GetFormalCharge(m)
    raise SystemExit(f"{path}: pteridine C2/N1/N3 not found")

new, old = tautomer(sys.argv[1]), tautomer(sys.argv[2])
print(f"   (H on N1, H on N3, net charge): FOL-3H {new}, v4 FOL {old}")
sys.exit(0 if new == (0, 1, -2) and old == (1, 0, -2) else 1)
EOF
echo "== package FOL-3H ready  $(date)"
