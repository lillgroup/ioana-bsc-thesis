#!/bin/bash
#SBATCH --job-name=xw-build
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --array=0-6
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/build-%A_%a.out
# Complex REST2 rebuilt WITH the deposited crystal waters (md-tools' protocol: they are kept, moved
# after the solute, and the box is filled around them). Everything else is the run it replaces,
# copied unchanged: package (reused, never re-parameterised), build-top.config, scaler.config,
# REST2.config, protein + ligand coordinates (scripts/prep.py writes the same records as the old prep.py, plus the waters).
#   0 FOL-4LRH  <- FOL-REST2-v4/4LRH            3 4ibk <- REST2-conformers-v4/4ibk-complex
#   1 FOL-3QL3  <- FOL-REST2-v4/3QL3            4 4j7d <- REST2-conformers-v4/4j7d-complex
#   2 4ibi      <- REST2-conformers-v4/4ibi-complex  5 5mwh <- REST2-conformers-complex/5mwh
#                                                    6 5ov8 <- REST2-conformers-complex/5ov8
# Writes REST2-run1 (production, as in REST2.config) and test/REST2-test (300 exchanges).
set -uo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
NAMES=(FOL-4LRH FOL-3QL3 4ibi 4ibk 4j7d 5mwh 5ov8)
PKGS=(FOL FOL 4ibi 4ibk 4j7d 5mwh-water 5ov8-water)
LIGS=(FOL FOL 1D8 1DE I31 UWX AXN)
i=$SLURM_ARRAY_TASK_ID; n=${NAMES[$i]}; pkg=$ROOT/packages/${PKGS[$i]}; S=$ROOT/$n
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$(/scicore/home/lill0000/popa0002/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
echo "== $n  $(md-openmm --version)  $(git -C /scicore/home/lill0000/popa0002/software/MD-tools-protonation log --oneline -1)  on $(hostname)  $(date)"
fail() { echo "FAILED: $1"; exit 1; }

echo "== 1. prep: chain A + ligand(s) + crystal waters"
python -I "$ROOT/scripts/prep.py" "$ROOT" "$n" || fail prep

echo "== 2. build-top"
b=$S/build
( cd "$b" && md-openmm build-top --overwrite -i ../prep/complex.pdb --config build-top.config \
    -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) || { grep -v '^$' "$b/built.out" | tail -25; fail build-top; }
grep -E 'missing atoms|reused|created|retained|net charge|atoms  |waters  |solvation|status' "$b/built.out" | head -14
# complex builds report the mapped package as "resname 'X' -> <compound>/<param_id>"
PID=$(awk '/^parameter_id:/{print $2}' "$pkg/parameter.config")
grep -q "resname '${LIGS[$i]}' -> .*/$PID" "$b/built.out" || fail "ligand not mapped to package $PID"
echo "   ${LIGS[$i]} -> $PID (the committed package, reused)"

echo "== 3. ligand residue = the scaler's mask"
IDX=$(grep -oE 'mask: ":[0-9]+"' "$b/scaler.config" | head -1 | grep -oE '[0-9]+')
RNAME=$(python -I - "$b/built.pdb" "$IDX" <<'EOF'
import sys
from openmm.app import PDBFile
top = PDBFile(sys.argv[1]).topology
r = list(top.residues())[int(sys.argv[2]) - 1]
hoh = sum(x.name == "HOH" for x in top.residues())
print(r.name)
print(f"   residue {sys.argv[2]} = {r.name}; {hoh} waters in the box", file=sys.stderr)
EOF
) || fail "residue check"
[ "$RNAME" = "${LIGS[$i]}" ] || [ "$RNAME" = "$(basename "$(ls "$pkg"/L0?.pdb 2>/dev/null || echo none)" .pdb)" ] \
  || fail "residue :$IDX is $RNAME, expected ${LIGS[$i]}"
cp -p "$pkg/molecule.sdf" "$b/$RNAME.sdf"     # bond orders for the scaler

echo "== 4. scaler"
( cd "$b" && md-openmm build-top --overwrite --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1 ) \
  || { tail -20 "$b/scaler.out"; fail scaler; }
grep -E 'selection|unscaled torsions|amide' "$b/scaler.out" | head -5

echo "== 5. build-md: production and a 300-exchange test"
( cd "$S" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
  || { tail -20 "$S/build-md.out"; fail build-md; }
mkdir -p "$S/test"; ln -sfn ../build "$S/test/build"
sed -E 's/^(  number_of_exchanges:) *[0-9]+.*/\1 300   # TEST: 0.6 ns per state/' "$S/REST2.config" > "$S/test/REST2.config"
( cd "$S/test" && md-openmm build-md -odir ./REST2-test --config REST2.config > build-md.out 2>&1 ) \
  || { tail -20 "$S/test/build-md.out"; fail "build-md (test)"; }
grep -h 'number_of_exchanges\|number_of_replicas' "$S/input/REST2.in" "$S/test/input/REST2.in"
echo "== $n built  $(date)"
