#!/bin/bash
#SBATCH --job-name=report-cut
#SBATCH --partition=scicore
#SBATCH --qos=1day
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/report-%x-%j.out
# The clustering figures of ONE compound with Shu-Yu's TORSION PARTITIONING (barrier cuts,
# ../shuyu_clustering/cluster_cut.py, MD-tools main 306b109, env analysis-env-cut), water and every
# complex, on the same inputs as report.sbatch: prepare.py's trajectories and torsion files, so t<k> is
# the same bond as in the symmetry-first figures. MD-tools' defaults, every frame (stride 1).
# Figures -> results/<compound>/<version>/<compound>_<system>_clusters_cut.png (+ _cut_split_tree.png,
# _cut_free_energy.png, _cut_groups.csv); a section is added to that folder's README.
#   sbatch -J report-cut-4ibi analysis/report/report_cut.sbatch 4ibi
set -uo pipefail
NAME=${1:?compound name from compounds.csv}
T=/scicore/home/lill0000/popa0002/thesis
REPO=$T/bsc-thesis
HERE=$REPO/analysis/report
WORK=$T/analysis-report/$NAME
M=$WORK/manifest.json
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
fail() { echo "FAILED: $1"; exit 1; }
if [ ! -s "$M" ]; then                       # inputs not prepared yet: same step as report.sbatch
  micromamba activate analysis-env-main
  python -I "$HERE/prepare.py" --name "$NAME" || fail prepare
fi
micromamba activate analysis-env-cut
echo "== report-cut $NAME  md-tools $(python -c 'import md_tools; print(md_tools.__version__)') ($(git -C /scicore/home/lill0000/popa0002/software/MD-tools-cut log --oneline -1 | cut -c1-7))  on $(hostname)  $(date)"
VERSION=$(python -I -c "import json; print(json.load(open('$M'))['version'])")
OUTDIR=$REPO/results/$NAME/$VERSION
mkdir -p "$OUTDIR"
NC=$(python -I -c "import json; print(len(json.load(open('$M'))['complexes']))")
summary=()
for KEY in water $(seq -s ' ' -f 'complexes.%g' 0 $((NC - 1)) 2>/dev/null); do
  [ "$KEY" = "complexes.-1" ] && continue
  mapfile -t A < <(python -I - "$M" "$KEY" "$NAME" <<'EOF'
import json, os, sys
m = json.load(open(sys.argv[1])); key = sys.argv[2]
s = m["water"] if key == "water" else m["complexes"][int(key.split(".")[1])]
print(s["label"])
args = ["--traj", s["traj"], "--top", s["top"], "--sdf", s["sdf"],
        "--name", f"{sys.argv[3]} {'water' if key == 'water' else s['label'] + ' complex'}"]
for k in ("torsions", "torsions_rigid", "torsions_polarH"):
    if os.path.getsize(f"{s['torsions']}/{k}.txt") > 0:
        args += ["--torsions-file", f"{s['torsions']}/{k}.txt"]
for r in s["refs"]:
    args += ["--ref-sdf", r]
if m["max_matches"]:
    args += ["--max-matches", str(m["max_matches"])]
if s.get("rest2_config") and os.path.exists(s["rest2_config"]):
    args += ["--rest2-config", s["rest2_config"]]
print("\n".join(args))
EOF
)
  L=${A[0]}; OUT=$WORK/$L-cut
  python "$REPO/analysis/shuyu_clustering/cluster_cut.py" "${A[@]:1}" --out "$OUT" > "$OUT.log" 2>&1 \
    || { tail -20 "$OUT.log"; fail "cluster_cut $L"; }
  line=$(grep -E '^[0-9]+ clusters; masses' "$OUT.log")
  short=$(python -I -c "import ast,sys; m=ast.literal_eval(sys.argv[1].split(': ',1)[1]); print(', '.join(map(str, m[:10])) + (f' ... (+{len(m) - 10} more, {sum(m[10:]):.1f} % together)' if len(m) > 10 else ''))" "$line")
  echo "   $L: ${line%%;*}: $short"; summary+=("| $L | ${line%%;*} | $short |")
  cp "$OUT/clusters_overview.png" "$OUTDIR/${NAME}_${L}_clusters_cut.png"
  cp "$OUT/groups.csv" "$OUTDIR/${NAME}_${L}_cut_groups.csv"
  [ -s "$OUT/split_tree.png" ] && cp "$OUT/split_tree.png" "$OUTDIR/${NAME}_${L}_cut_split_tree.png"
  cp "$OUT/clusters_free_energy.png" "$OUTDIR/${NAME}_${L}_cut_free_energy.png"
done
README=$OUTDIR/README.md
[ -s "$README" ] && python -I - "$README" <<'EOF'
import sys
p = sys.argv[1]; s = open(p).read()
i = s.find("\n## Torsion partitioning")
open(p, "w").write(s if i < 0 else s[:i + 1])
EOF
{
  echo; echo "## Torsion partitioning (barrier cuts)"; echo
  echo "\`*_clusters_cut.png\` (+ \`_cut_split_tree.png\`, \`_cut_free_energy.png\`, \`_cut_groups.csv\`): Shu-Yu's"
  echo "\`cut_clustering\` (MD-tools main $(git -C /scicore/home/lill0000/popa0002/software/MD-tools-cut log --oneline -1 | cut -c1-7)), \`analysis/report/report_cut.sbatch\` (Slurm job ${SLURM_JOB_ID:-}, $(date +%F));"
  echo "MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above."
  echo; echo "| system | clusters | masses (%) |"; echo "|---|---|---|"
  printf '%s\n' "${summary[@]}"
} >> "$README"
echo "== $NAME done -> $OUTDIR  $(date)"
