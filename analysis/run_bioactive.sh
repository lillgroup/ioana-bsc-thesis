#!/usr/bin/env bash
# Bioactive-conformer analysis (shuyu_clustering/bioactive_rmsd.py): RMSD of the ligand conformation to
# its crystal for the whole water ensemble (RMSD-wat) and per cluster (RMSD-wat-A, ...), and the same for
# the protein-bound ensemble (RMSD-prot, RMSD-prot-A, ...). v4 runs, longest run each, state 0.
# Clusterings (labels.csv) used:
#   folate         water: symfirst-0651-gauge (--gauge torsions_3); 4LRH / 3QL3 complexes: the same
#                  setting (run_symfirst_0651.sh with TAGBASE=symfirst-0651-gauge FIT_ARGS="--gauge torsions_3")
#   4ibi 4ibk 4j7d symfirst-0651-rec (recursive split had no effect on these: = the plain clustering);
#                  4ibi complex: -rec-single (one state, the allow_single_cluster=True fallback)
# Water ensembles are compared with the chain-A crystal of the complex. Outputs in results/_incoming/bioactive/ (then sorted into results/<compound>/<version>/bioactive/).
# Written for the local WSL machine; on sciCORE submit the python calls with sbatch.
#   bash run_bioactive.sh            # everything
#   bash run_bioactive.sh 4j7d       # one system (folate-4LRH folate-3QL3 4ibi 4ibk 4j7d)
set -uo pipefail
T=${THESIS:-$HOME/thesis-sbrinz-data/thesis}
REPO=${REPO:-$HOME/bsc-thesis}
OUT=$REPO/results/_incoming/bioactive
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/software/md-stack}
PY=("$MAMBA_ROOT_PREFIX/bin/micromamba" run -n analysis-env-main python -I)
mkdir -p "$OUT"
ONLY=${1:-}
want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }
run() {   # run <out name> [bioactive_rmsd.py args...]
  local N=$1; shift
  echo "== $N  $(date +%T)"
  "${PY[@]}" "$REPO/analysis/shuyu_clustering/bioactive_rmsd.py" "$@" --out "$OUT/$N" 2>&1 | grep -iv warning > "$OUT/$N.log"
  grep -E '^\| ' "$OUT/$N.log"
}

F=$T/FOL-REST2-v4
for S in 4LRH 3QL3; do
  want "folate-$S" || continue
  run "folic-acid-v4_$S" --name "Folic acid (v4)" --sdf "$F/water/build/FOL.sdf" --top "$F/water/build/built.solute.pdb" \
    --crystal "$REPO/crystal/fol_${S}_A.sdf" --crystal-name "$S" \
    --ensemble wat "$F/water/REST2-run1/solute_state0_prod1.nc" \
               "$F/water/folic-acid-v4-water-clustering-symfirst-0651-gauge/labels.csv" \
    --ensemble prot "$F/$S/analysis-shuyu/FOL_state0_prod1.nc" \
               "$F/$S/analysis-shuyu/folic-acid-v4-$S-complex-clustering-symfirst-0651-gauge/labels.csv"
done

C=$T/REST2-conformers-v4
for spec in 4ibi:L01 4ibk:L02 4j7d:L03; do
  IFS=: read -r n res <<< "$spec"
  want "$n" || continue
  W=$C/$n-water; X=$C/$n-complex/analysis-symfirst
  WL=$W/analysis-symfirst/$n-v4-water-clustering-symfirst-0651-rec/labels.csv
  XL=$X/$n-v4-complex-clustering-symfirst-0651-rec/labels.csv
  [ -s "$XL" ] || XL=$X/$n-v4-complex-clustering-symfirst-0651-rec-single/labels.csv
  for crys in "$REPO/crystal/${n}"_*_A.sdf; do
    run "$n-v4" --name "$n ligand (v4)" --sdf "$W/build/$res.sdf" --top "$W/build/built.solute.pdb" \
      --crystal "$crys" --crystal-name "$n chain A" \
      --ensemble wat "$W/REST2-run1/solute_state0_prod1.nc" "$WL" \
      --ensemble prot "$X/${res}_state0_prod1.nc" "$XL"
  done
done
echo "== all done  $(date +%T); results in $OUT"
