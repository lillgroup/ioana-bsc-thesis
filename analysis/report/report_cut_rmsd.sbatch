#!/bin/bash
#SBATCH --job-name=report-cut-rmsd
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/report-%x-%j.out
# Inter-cluster RMSD + cluster network for the BARRIER-CUT clusters of ONE compound (water and every
# complex), on the finished report_cut.sbatch runs ($WORK/<system>-cut): ../shuyu_clustering/
# cluster_cut_rmsd.py, clusters >= 1 % of the frames. Nothing is re-clustered.
# Figures -> results/<compound>/<version>/<compound>_<system>_cut_intercluster_rmsd.png, _cut_cluster_network.png
#   sbatch -J report-cut-rmsd-4ibi analysis/report/report_cut_rmsd.sbatch 4ibi
set -uo pipefail
NAME=${1:?compound name from compounds.csv}
T=/scicore/home/lill0000/popa0002/thesis
REPO=$T/bsc-thesis
WORK=$T/analysis-report/$NAME
M=$WORK/manifest.json
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-cut
export OMP_NUM_THREADS=1
fail() { echo "FAILED: $1"; exit 1; }
echo "== report-cut-rmsd $NAME  on $(hostname)  $(date)"
VERSION=$(python -I -c "import json; print(json.load(open('$M'))['version'])")
OUTDIR=$REPO/results/$NAME/$VERSION
while IFS='|' read -r L TRAJ TOP SDF MM SYS; do
  OUT=$WORK/$L-cut
  [ -s "$OUT/labels.csv" ] || fail "no barrier-cut run in $OUT (report_cut.sbatch first)"
  python "$REPO/analysis/shuyu_clustering/cluster_cut_rmsd.py" --cut-dir "$OUT" --traj "$TRAJ" --top "$TOP" \
    --sdf "$SDF" ${MM:+--max-matches $MM} --rmsd-workers "${SLURM_CPUS_PER_TASK:-16}" --name "$NAME $SYS" \
    > "$OUT-rmsd.log" 2>&1 || { tail -20 "$OUT-rmsd.log"; fail "cluster_cut_rmsd $L"; }
  head -1 "$OUT-rmsd.log" | sed 's/^/   /'
  cp "$OUT/intercluster_rmsd.png" "$OUTDIR/${NAME}_${L}_cut_intercluster_rmsd.png"
  cp "$OUT/cluster_network.png" "$OUTDIR/${NAME}_${L}_cut_cluster_network.png"
done < <(python -I - "$M" <<'EOF'
import json, sys
m = json.load(open(sys.argv[1]))
for s in [m["water"]] + m["complexes"]:
    sys_ = "water" if s is m["water"] else s["label"] + " complex"
    print("|".join([s["label"], s["traj"], s["top"], s["sdf"], str(m["max_matches"] or ""), sys_]))
EOF
)
README=$OUTDIR/README.md
grep -q "_cut_intercluster_rmsd.png" "$README" 2>/dev/null || cat >> "$README" <<EOF

\`*_cut_intercluster_rmsd.png\`, \`*_cut_cluster_network.png\`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
\`analysis/report/report_cut_rmsd.sbatch\` (\`shuyu_clustering/cluster_cut_rmsd.py\`).
EOF
echo "== $NAME done -> $OUTDIR  $(date)"
