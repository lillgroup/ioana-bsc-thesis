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
# Optional 2nd argument: the lowest barrier cut kT_cut (E*, MD-tools default 1). Any other value runs
# into <system>-cut-kT<X>/ and names the figures *_clusters_cut_kT<X>.png etc., next to the default ones:
#   sbatch -J report-cut-kT2-4ibi analysis/report/report_cut.sbatch 4ibi 2
set -uo pipefail
NAME=${1:?compound name from compounds.csv}
KT=${2:-1}
if [ "$KT" = 1 ]; then SUF=""; KTARG=(); else SUF="_kT$KT"; KTARG=(--kT-cut "$KT"); fi
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
  L=${A[0]}; OUT=$WORK/$L-cut${SUF/_/-}
  python "$REPO/analysis/shuyu_clustering/cluster_cut.py" "${A[@]:1}" "${KTARG[@]}" --out "$OUT" > "$OUT.log" 2>&1 \
    || { tail -20 "$OUT.log"; fail "cluster_cut $L"; }
  line=$(grep -E '^[0-9]+ clusters; masses' "$OUT.log")
  short=$(python -I -c "import ast,sys; m=ast.literal_eval(sys.argv[1].split(': ',1)[1]); print(', '.join(map(str, m[:10])) + (f' ... (+{len(m) - 10} more, {sum(m[10:]):.1f} % together)' if len(m) > 10 else ''))" "$line")
  echo "   $L: ${line%%;*}: $short"; summary+=("| $L | ${line%%;*} | $short |")
  cp "$OUT/clusters_overview.png" "$OUTDIR/${NAME}_${L}_clusters_cut$SUF.png"
  cp "$OUT/groups.csv" "$OUTDIR/${NAME}_${L}_cut${SUF}_groups.csv"
  [ -s "$OUT/split_tree.png" ] && cp "$OUT/split_tree.png" "$OUTDIR/${NAME}_${L}_cut${SUF}_split_tree.png"
  cp "$OUT/clusters_free_energy.png" "$OUTDIR/${NAME}_${L}_cut${SUF}_free_energy.png"
done
README=$OUTDIR/README.md
if [ "$KT" = 1 ]; then HEAD="## Torsion partitioning (barrier cuts)"; else HEAD="## Torsion partitioning (barrier cuts, kT_cut $KT)"; fi
[ -s "$README" ] && python -I - "$README" "$HEAD" <<'EOF'
import sys
p, head = sys.argv[1], sys.argv[2]; s = open(p).read()
i = s.find("\n" + head + "\n")                     # this kT's section only, up to the next section
if i >= 0:
    j = s.find("\n## ", i + 1)
    s = s[:i] + (s[j:] if j >= 0 else "\n")
open(p, "w").write(s)
EOF
{
  echo; echo "$HEAD"; echo
  echo "\`*_clusters_cut$SUF.png\` (+ \`_cut${SUF}_split_tree.png\`, \`_cut${SUF}_free_energy.png\`, \`_cut${SUF}_groups.csv\`): Shu-Yu's"
  echo "\`cut_clustering\` (MD-tools main $(git -C /scicore/home/lill0000/popa0002/software/MD-tools-cut log --oneline -1 | cut -c1-7)), \`analysis/report/report_cut.sbatch\` (Slurm job ${SLURM_JOB_ID:-}, $(date +%F));"
  echo "MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut $KT (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above."
  echo; echo "| system | clusters | masses (%) |"; echo "|---|---|---|"
  printf '%s\n' "${summary[@]}"
} >> "$README"
echo "== $NAME done -> $OUTDIR  $(date)"
