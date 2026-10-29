#!/usr/bin/env bash
# Shu-Yu's symmetry-first clustering with MD-tools 0.6.5.1 (main 2e9e205: torsional folding of
# resonance-equivalent atoms, crystal -> cluster assignment, symmetry-corrected inter-cluster RMSD),
# shuyu_clustering/cluster_prune.py + cluster_convergence.py, on the v4 runs only, longest run each:
#   paracetamol-v4 (10 ns)   folate water (100 ns) / 4LRH / 3QL3 (50 ns)
#   4ibi 4ibk 4j7d, water (60 ns) and complex (40 ns)
# Same inputs as run_symfirst.sbatch / run_small_symfirst.sbatch (state 0, torsion files, stride,
# windows, crystal refs). Setting, by the user's decision (2026-10-26):
#   TAG symfirst-0651         allow_single_cluster=False (MD-tools' default since 0.6.5.1), always
#   TAG symfirst-0651-single  allow_single_cluster=True, ONLY as a fallback when the default finds
#                             no cluster at all (NO_CLUSTER.txt: a one-state ensemble, e.g. the 4ibi
#                             complex) -- never run alongside a default result
# Outputs next to each run (<name>-clustering-<TAG>/, <name>-convergence-<TAG>/, logs); the figures
# are collected in results/symfirst-0651/ (incl. cluster_network.png, shuyu_clustering/cluster_network.py). A step whose output already exists is skipped.
# Written for the local WSL machine (no Slurm, 7 GB RAM: one job at a time). On a cluster, wrap
# each python call in srun / submit it with sbatch instead.
#   bash run_symfirst_0651.sh                 # everything
#   bash run_symfirst_0651.sh 4ibi-v4-water   # one system (names below)
#   STEPS=both bash run_symfirst_0651.sh      # also the convergence check (default: clustering only,
#                                             # by the user's decision; STEPS=conv: convergence only)
#   FORCE=1 ...                               # redo steps whose output exists (e.g. after a figure change)
#   TAGBASE=symfirst-0651-rec FIT_ARGS="--mass-floor 0.005 --k 20 --min-vote-fraction 0.95 --recursive-max-rmsd 2.0" \
#     bash run_symfirst_0651.sh 4ibi-v4-water   # a variant, in its own folders and results/<TAGBASE>/
# Convergence re-clusters every cumulative 5 ns window (10-20 fits per system), so it costs 10-20x a
# clustering; running all clusterings first gets the results out before the slow part.
set -uo pipefail
T=${THESIS:-$HOME/thesis-sbrinz-data/thesis}
STEPS=${STEPS:-cluster}
REPO=${REPO:-$HOME/bsc-thesis}
SHUYU=$REPO/shuyu_clustering
CRYSTAL=$REPO/crystal
TAGBASE=${TAGBASE:-symfirst-0651}        # output tag; e.g. symfirst-0651-rec for a variant
read -r -a FIT_ARGS <<< "${FIT_ARGS:-}"   # extra cluster_prune.py settings, e.g. "--recursive-max-rmsd 2.0"
COLLECT=$REPO/results/$TAGBASE
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/software/md-stack}
PY=("$MAMBA_ROOT_PREFIX/bin/micromamba" run -n analysis-env-main python)
mkdir -p "$COLLECT"

run_system() {   # name root traj step cfg [args...] -- uses ARGS / REFS / EXTRA set by the caller
  local NAME=$1 ROOT=$2 TRAJ=$3 STEP=$4 CFG=$5
  cd "$ROOT" || return 1
  local TAGS=("$TAGBASE") k=0
  while [ $k -lt ${#TAGS[@]} ]; do
    local TAG=${TAGS[$k]}; k=$((k + 1))
    local FLAG=(); [ "$TAG" = "$TAGBASE-single" ] && FLAG=(--allow-single-cluster)
    local COMMON=(--traj "$TRAJ" "${ARGS[@]}" --name "$NAME" --rest2-config "$CFG" "${FLAG[@]}" "${FIT_ARGS[@]}")
    local OUT=$NAME-clustering-$TAG
    if [ "$STEPS" != conv ]; then
    echo "== $NAME $TAG clustering  $(date +%T)"
    if [ -z "${FORCE:-}" ] && { [ -s "$OUT/summary.json" ] || [ -s "$OUT/NO_CLUSTER.txt" ]; }; then
      echo "   exists, skipped"
    else
      "${PY[@]}" "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" "${EXTRA[@]}" --out "$OUT" > "$OUT.log" 2>&1
      rc=$?
      [ $rc = 3 ] && echo "   NO CLUSTER (one state at allow_single_cluster=False)"
      [ $rc != 0 ] && [ $rc != 3 ] && { echo "   FAILED ($rc)"; tail -5 "$OUT.log"; }
    fi
    grep -hE '^n_clusters|^masses|symmetry operations acting|^crystal .*best cluster' "$OUT.log" | sed 's/^/   /'
    [ -s "$OUT/intercluster_rmsd.json" ] && "${PY[@]}" "$SHUYU/cluster_network.py" "$OUT" --name "$NAME ($TAG)" > /dev/null 2>&1
    for f in clusters_overview.png intercluster_rmsd.png crystal_assignment.csv intercluster_rmsd.csv cluster_network.png; do
      [ -s "$OUT/$f" ] && cp "$OUT/$f" "$COLLECT/${NAME}_${TAG}_$f"
    done
    fi
    # the default found no cluster (one state): only then fall back to allow_single_cluster=True
    if [ "$TAG" = "$TAGBASE" ] && [ -s "$OUT/NO_CLUSTER.txt" ]; then
      echo "   -> one state: falling back to allow_single_cluster=True"
      TAGS+=("$TAGBASE-single")
    fi
    [ "$STEPS" = cluster ] && continue
    local CONV=$NAME-convergence-$TAG
    echo "== $NAME $TAG convergence  $(date +%T)"
    if [ -s "$CONV/convergence.csv" ]; then
      echo "   exists, skipped"
    else
      "${PY[@]}" "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns "$STEP" --plateau 3 --out "$CONV" > "$CONV.log" 2>&1 \
        || { echo "   FAILED"; tail -5 "$CONV.log"; }
    fi
    grep -h converged "$CONV.log" | tail -1 | sed 's/^/   /'
    [ -s "$CONV/convergence.png" ] && cp "$CONV/convergence.png" "$COLLECT/${NAME}_${TAG}_convergence.png"
  done
}

fol_args() {   # the water run's build + torsion files (complexes: same 51 atoms, same order)
  ARGS=(--sdf "$T/FOL-REST2-v4/water/build/FOL.sdf" --top "$T/FOL-REST2-v4/water/build/built.solute.pdb"
        --torsions-file "$T/FOL-REST2-v4/water/REST2-run1/torsions.txt"
        --torsions-file "$T/FOL-REST2-v4/water/REST2-run1/torsions_rigid.txt"
        --torsions-file "$T/FOL-REST2-v4/water/REST2-run1/torsions_polarH.txt" --stride 5)
  EXTRA=(--pair-2d 5,3 --pair-2d 4,3 --pair-2d 2,3)
}
small_args() {   # small_args <name> <res> <maxmatch> water|complex
  local n=$1 res=$2 mm=$3 kind=$4 W=$T/REST2-conformers-v4/$1-water
  ARGS=(--top "$W/build/built.solute.pdb" --sdf "$W/build/$res.sdf"
        --torsions-file "$W/REST2-run1/torsions.txt" --torsions-file "$W/REST2-run1/torsions_rigid.txt"
        --torsions-file "$W/REST2-run1/torsions_polarH.txt" --stride 5 --max-matches "$mm")
  REFS=(); EXTRA=()
  if [ "$kind" = water ]; then for f in "$CRYSTAL/${n}"_*.sdf; do REFS+=(--ref-sdf "$f"); done
  else for f in "$CRYSTAL/${n}"_*_A.sdf; do REFS+=(--ref-sdf "$f"); done; fi
}

ONLY=${1:-}
want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }

# small molecules first (fast), folate last (10 000-frame dense matrices)
for spec in 4ibi:L01:200000 4ibk:L02:200000 4j7d:L03:8000000; do
  IFS=: read -r n res mm <<< "$spec"
  if want "$n-v4-water"; then
    small_args "$n" "$res" "$mm" water
    mkdir -p "$T/REST2-conformers-v4/$n-water/analysis-symfirst"
    run_system "$n-v4-water" "$T/REST2-conformers-v4/$n-water/analysis-symfirst" \
      "$T/REST2-conformers-v4/$n-water/REST2-run1/solute_state0_prod1.nc" 5 "$T/REST2-conformers-v4/$n-water/REST2.config"
  fi
  if want "$n-v4-complex"; then
    small_args "$n" "$res" "$mm" complex
    run_system "$n-v4-complex" "$T/REST2-conformers-v4/$n-complex/analysis-symfirst" \
      "${res}_state0_prod1.nc" 5 "$T/REST2-conformers-v4/$n-complex/REST2.config"
  fi
done
if want paracetamol-v4; then
  ARGS=(--sdf build/TYL.sdf --top build/built.solute.pdb --stride 1); REFS=(); EXTRA=()
  run_system paracetamol-v4 "$T/PARA-REST2-v4" REST2-run1/solute_state0_prod1.nc 0.5 REST2.config
fi
if want folic-acid-v4-water; then
  fol_args; REFS=(--ref-sdf "$CRYSTAL/fol_4LRH_A.sdf" --ref-sdf "$CRYSTAL/fol_3QL3_A.sdf")
  run_system folic-acid-v4-water "$T/FOL-REST2-v4/water" REST2-run1/solute_state0_prod1.nc 5 REST2.config
fi
for S in 4LRH 3QL3; do
  if want "folic-acid-v4-$S-complex"; then
    fol_args; REFS=(--ref-sdf "$CRYSTAL/fol_${S}_A.sdf")
    run_system "folic-acid-v4-$S-complex" "$T/FOL-REST2-v4/$S/analysis-shuyu" FOL_state0_prod1.nc 5 ../REST2.config
  fi
done
echo "== all done  $(date +%T); figures in $COLLECT"
