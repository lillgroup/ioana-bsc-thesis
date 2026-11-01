#!/bin/bash
#SBATCH --job-name=rb-build
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-14
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds/logs/build-%A_%a.out
# REST2-rotbonds: one compound per array task (row of compounds.csv, 0-based), CPU only.
#   1. parameterise the ligand (new compound -> new package, the v4 route of REST2-conformers-v4:
#      neutral input, Open Babel pH 7.4, Sage + AM1-BCC); skipped if packages/<n> exists
#   2. water build (as REST2-conformers-v4/<n>-water): REST2-run1 (60 ns) + test/REST2-test (300 exchanges)
#   3. complex: prep.py (chain + ligand + crystal waters), build-top (package reused), the ligand
#      residue number into scaler.config / REST2.config, scaler, REST2-run1 (40 ns) + test
# A part whose REST2-run1/run.sh exists is not rebuilt (never overwrites a started run).
#   sbatch scripts/build.sbatch                 all 15        sbatch --array=3 scripts/build.sbatch   only 4p6w
set -uo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds
i=$SLURM_ARRAY_TASK_ID
IFS=, read -r n pdb lig chain res _ < <(sed -n "$((i + 2))p" "$ROOT/compounds.csv")
[ -n "$n" ] || { echo "no row $i in compounds.csv"; exit 2; }
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
echo "== $i $n ($pdb chain $chain, $lig -> package residue $res)  $(md-openmm --version)" \
     "$(git -C /scicore/home/lill0000/popa0002/software/MD-tools-protonation log --oneline -1)  on $(hostname)  $(date)"
fail() { echo "FAILED: $1"; exit 1; }
test_run() {   # test_run <system dir>: a 300-exchange copy of REST2.config sharing build/
  local S=$1
  mkdir -p "$S/test"; ln -sfn ../build "$S/test/build"
  sed -E 's/^(  number_of_exchanges:) *[0-9]+.*/\1 300   # TEST: 0.6 ns per state/' "$S/REST2.config" > "$S/test/REST2.config"
  ( cd "$S/test" && md-openmm build-md -odir ./REST2-test --config REST2.config > build-md.out 2>&1 ) \
    || { tail -20 "$S/test/build-md.out"; fail "build-md (test) $S"; }
}

echo "== 1. parameterise at pH 7.4"
cd "$ROOT/packages"
if [ -f "$n/parameter.config" ]; then
  echo "   package exists, reused"
else
  md-openmm build-top --parameterize -i "$n.input.sdf" --resname "$res" --config "$n.parameterize.config" \
      -op "$n/$res.pdb" -os "$n/$res.xml" -log "$n.parameterize.log" > "$n.parameterize.out" 2>&1 \
    || { grep -v '^$' "$n.parameterize.out" | tail -20; fail parameterise; }
fi
grep -E 'net_formal_charge|canonical_smiles|parameter_id' "$n/parameter.config"
grep -iE 'protonation:' "$n.parameterize.out" 2>/dev/null | head -3
PID=$(awk '/^parameter_id:/{print $2}' "$n/parameter.config")

echo "== 2. water"
S=$ROOT/$n-water; b=$S/build
if [ -x "$S/REST2-run1/run.sh" ]; then
  echo "   already built (REST2-run1/run.sh exists)"
else
  ( cd "$b" && md-openmm build-top --overwrite -i "../../packages/$n/$res.sdf" --config build-top.config \
      -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) || { tail -20 "$b/built.out"; fail "water build-top"; }
  grep -E 'reused|net charge|atoms  |protonation|status' "$b/built.out" | head -8
  cp -p "$ROOT/packages/$n/molecule.sdf" "$b/$res.sdf"
  ( cd "$b" && md-openmm build-top --overwrite --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1 ) \
    || { tail -20 "$b/scaler.out"; fail "water scaler"; }
  grep -E 'selection|unscaled torsions|amide' "$b/scaler.out" | head -5
  ( cd "$S" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
    || { tail -20 "$S/build-md.out"; fail "water build-md"; }
  test_run "$S"
fi
grep -hE 'number_of_replicas|number_of_exchanges' "$S/input/REST2.in" "$S/test/input/REST2.in"

echo "== 3. complex"
S=$ROOT/$n-complex; b=$S/build
if [ -x "$S/REST2-run1/run.sh" ]; then
  echo "   already built (REST2-run1/run.sh exists)"
else
  python -I "$ROOT/scripts/prep.py" "$ROOT" "$n" || fail prep
  ( cd "$b" && md-openmm build-top --overwrite -i ../prep/complex.pdb --config build-top.config \
      -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) || { grep -v '^$' "$b/built.out" | tail -25; fail "complex build-top"; }
  grep -E 'missing atoms|reused|created|retained|net charge|atoms  |waters  |solvation|WARNING|status' "$b/built.out" | head -14
  grep -q "resname '$lig' -> .*/$PID" "$b/built.out" || fail "ligand $lig not mapped to package $PID"
  echo "   $lig -> $PID (this compound's package, as in the water run)"
  read -r IDX RNAME < <(python -I - "$b/built.pdb" "$lig" "$res" <<'EOF'
import sys
from openmm.app import PDBFile
top = PDBFile(sys.argv[1]).topology
hits = [(k, r.name) for k, r in enumerate(top.residues(), 1) if r.name in (sys.argv[2], sys.argv[3])]
assert len(hits) == 1, hits
print(*hits[0])
print(f"   ligand residue :{hits[0][0]} ({hits[0][1]}); "
      f"{sum(r.name == 'HOH' for r in top.residues())} waters in the box", file=sys.stderr)
EOF
  ) || fail "ligand residue not found"
  [ -n "${IDX:-}" ] || fail "ligand residue not found"
  sed -i -E "s/mask: \":[A-Z0-9]+\"/mask: \":$IDX\"/; s/residue :[A-Z0-9]+ =/residue :$IDX =/" "$b/scaler.config" "$S/REST2.config"
  cp -p "$ROOT/packages/$n/molecule.sdf" "$b/$RNAME.sdf"     # bond orders for the scaler
  ( cd "$b" && md-openmm build-top --overwrite --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1 ) \
    || { tail -20 "$b/scaler.out"; fail "complex scaler"; }
  grep -E 'selection|unscaled torsions|amide' "$b/scaler.out" | head -5
  ( cd "$S" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
    || { tail -20 "$S/build-md.out"; fail "complex build-md"; }
  test_run "$S"
fi
grep -hE 'number_of_replicas|number_of_exchanges' "$S/input/REST2.in" "$S/test/input/REST2.in"
grep -h 'mask' "$b/scaler.config"
echo "== $n built  $(date)"
