#!/usr/bin/env bash
# Collect ONLY the simulation outputs of the local water REST2 runs that the analysis scripts read,
# into one folder to copy to sciCORE (2026-11-02). Same relative paths as under thesis/, so a script
# only needs its root changed (THESIS=<bundle>). Files are hard links: no extra disk space locally.
#   bash envs/water_runs_bundle.sh [<thesis root> [<bundle>]]
# Kept per run: REST2.config, input/REST2.in, build/ (solute + full topology, ligand SDF and package,
# configs), and from REST2-run*/: every state's solute trajectory, state 0's whole-system trajectory,
# the torsion files, the exchange records (exchange.csv, rem.log, remd_records/, runstate, restart,
# resolved.config). For 5mwh / 5ov8 also analysis-0-60ns/solute_state0_0-60ns.nc (state 0 stitched from
# REST2-run1 + REST2-run1-ext1) with its torsion files.
# Left out: clustering / convergence outputs, min/, eq/, checkpoints, rank/, run scripts, the
# replica-ordered REST2.nc / REST2.solute.nc (same frames as the state trajectories), states 1..n of
# the whole-system trajectory.
set -euo pipefail
SRC=$(cd "${1:-$HOME/thesis-sbrinz-data/thesis}" && pwd)
DEST=${2:-$HOME/water-runs-for-scicore}
RUNS=(
  FOL-REST2-v4/water FOL-REST2-mace/water
  REST2-conformers-v4/4ibi-water REST2-conformers-v4/4ibk-water REST2-conformers-v4/4j7d-water
  REST2-conformers/5mwh REST2-conformers/5ov8
  PARA-REST2-v4
  REST2-replicates/FOL/rep1 REST2-replicates/FOL/rep2 REST2-replicates/FOL/rep3 REST2-replicates/PARA/rep1
)
FILTER=(
  --include='/REST2.config' --include='/input/' --include='/input/REST2.in'
  --include='/build/' --include='/build/built.solute.pdb' --include='/build/built.pdb'
  --include='/build/*.sdf' --include='/build/*.config' --include='/build/ligands/***'
  --include='/build/parameter/***'
  --include='/REST2-run*/'
  --include='/REST2-run*/solute_state*_prod1.nc' --include='/REST2-run*/whole_state0_prod1.nc'
  --include='/REST2-run*/torsions*.txt' --include='/REST2-run*/exchange.csv'
  --include='/REST2-run*/rem.log' --include='/REST2-run*/remd_records/***'
  --include='/REST2-run*/REST2.runstate.json' --include='/REST2-run*/restart.json'
  --include='/REST2-run*/resolved.config'
  --include='/analysis-0-60ns/' --include='/analysis-0-60ns/solute_state0_0-60ns.nc'
  --include='/analysis-0-60ns/torsions*.txt'
  --exclude='*'
)
mkdir -p "$DEST"
for r in "${RUNS[@]}"; do
  [ -d "$SRC/$r" ] || { echo "missing $SRC/$r"; exit 1; }
  mkdir -p "$DEST/$r"
  # patterns are relative to the run folder; --link-dest makes the copies hard links
  rsync -a --prune-empty-dirs --link-dest="$SRC/$r" "${FILTER[@]}" "$SRC/$r/" "$DEST/$r/"
done
cat > "$DEST/README.txt" <<EOF
Water REST2 runs (simulation outputs only), collected from the local copy of the sbrinz data by
bsc-thesis/envs/water_runs_bundle.sh on $(date +%F). Paths as under thesis/ on sbrinz; see that
script for what was kept and left out, and bsc-thesis/SETUP.md section 5 for what each run is.
EOF
du -sh "$DEST"; find "$DEST" -type f | wc -l
