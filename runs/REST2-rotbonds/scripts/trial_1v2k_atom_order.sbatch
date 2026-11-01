#!/bin/bash
#SBATCH --job-name=rb-1v2k-order
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-2
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds/logs/trial-1v2k-order-%A_%a.out
# 1v2k cation: md-tools' parameterisation self-check refused the package because its improper
# torsions list the three outer atoms in the other (equivalent) cyclic order than openmmforcefields'
# reference. The neutral molecule, written in the same atom order, passed. Atom order in the input is
# arbitrary, so try the SAME molecule (same coordinates, same chemistry) in three other orders and
# see whether the two builds then agree. Trial packages go to packages/_trial/v<k>/, nothing else.
#   0 every H right after its heavy atom      1 RDKit canonical atom order      2 H block reversed
set -uo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds
k=$SLURM_ARRAY_TASK_ID
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
D=$ROOT/packages/_trial; mkdir -p "$D"
python -I - "$ROOT/packages/1v2k.input.sdf" "$D/1v2k-v$k.input.sdf" "$k" <<'EOF'
import sys
from rdkit import Chem
m = Chem.MolFromMolFile(sys.argv[1], removeHs=False); k = int(sys.argv[3])
heavy = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
hs = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 1]
if k == 0:
    order = []
    for i in heavy:
        order.append(i)
        order += [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1]
elif k == 1:
    ranks = list(Chem.CanonicalRankAtoms(m, breakTies=True))
    order = sorted(range(m.GetNumAtoms()), key=lambda i: ranks[i])
else:
    order = heavy + hs[::-1]
r = Chem.RenumberAtoms(m, order); r.SetProp("_Name", m.GetProp("_Name"))
assert Chem.MolToSmiles(r) == Chem.MolToSmiles(m) and Chem.GetFormalCharge(r) == 1
Chem.MolToMolFile(r, sys.argv[2])
print(f"variant {k}: same molecule, charge {Chem.GetFormalCharge(r)}, new order written")
EOF
cd "$D"
md-openmm build-top --parameterize -i "1v2k-v$k.input.sdf" --resname X02 --config ../1v2k.parameterize.config \
    -op "v$k/X02.pdb" -os "v$k/X02.xml" -log "v$k.parameterize.log" > "v$k.parameterize.out" 2>&1
rc=$?
grep -E 'Open Babel|PackageError' "v$k.parameterize.out" | cut -c1-150
[ $rc -eq 0 ] && { echo "VARIANT $k PASSED"; grep -E 'net_formal_charge|parameter_id' "v$k/parameter.config"; } || echo "VARIANT $k REFUSED"
