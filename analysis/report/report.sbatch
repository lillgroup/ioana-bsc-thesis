#!/bin/bash
#SBATCH --job-name=report
#SBATCH --partition=scicore
#SBATCH --qos=1day
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/report-%x-%j.out
# The figures of ONE compound (a row of compounds.csv), from its current water and complex runs:
#   1. prepare.py          water + complex ligand trajectories (state 0), one torsion numbering for all
#   2. cluster_prune.py    Shu-Yu's symmetry-first clustering (MD-tools 0.6.5.2), --no-single-cluster,
#                          water and every complex
#   3. rmsd_density.py     ligand conformation RMSD to the crystal, water vs complex
#   4. check_bound.py      complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q
#   5. water_vs_complex.py torsion distributions, water vs complex, one panel per torsion
# Figures -> bsc-thesis/results/<compound>/<version>/ (version from compounds.csv), with a README
# saying which runs and settings; working files in $THESIS/analysis-report/<compound>/.
# If a run of the compound has not finished, the job resubmits itself after that run's jobs.
#   sbatch -J report-4ibi analysis/report/report.sbatch 4ibi        (or: bash analysis/report/submit.sh)
set -uo pipefail
NAME=${1:?compound name from compounds.csv}
T=/scicore/home/lill0000/popa0002/thesis
REPO=$T/bsc-thesis
HERE=$REPO/analysis/report
SHUYU=$REPO/analysis/shuyu_clustering
WORK=$T/analysis-report/$NAME
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
fail() { echo "FAILED: $1"; exit 1; }
quiet() { "$@" 2>&1 | grep -iv warning; return "${PIPESTATUS[0]}"; }   # run, drop warnings, keep its status
echo "== report $NAME  md-tools $(python -c 'import importlib.metadata as m; print(m.version("md-tools"))')  on $(hostname)  $(date)"

# ---- 0. are the runs finished? ------------------------------------------------------------------
ROW=$(python -I -c "import csv,sys; r={x['name']:x for x in csv.DictReader(open('$HERE/compounds.csv'))}['$NAME']; print(r['water'], r['water_ext'] or '-', r['complexes'] or '-', r['version'])")
read -r WATER WEXT COMPLEXES VERSION <<< "$ROW"
RUNS=("$WATER/REST2-run1")
[ "$WEXT" != - ] && for e in ${WEXT//;/ }; do RUNS+=("$WATER/$e"); done
CRUNS=()
[ "$COMPLEXES" != - ] && for spec in ${COMPLEXES//;/ }; do IFS== read -r _ run _ <<< "$spec"; CRUNS+=("$run"); RUNS+=("$run/REST2-run1"); done
missing=()
for r in "${RUNS[@]}"; do
  f=$T/$r/remd_records/REST2_prod1.out
  [ -s "$T/$r/solute_state0_prod1.nc" ] && { [ ! -e "$f" ] || grep -q 'run_status: completed' "$f"; } || missing+=("$r")
done
if [ ${#missing[@]} -gt 0 ]; then
  jobs=$(squeue -u "$USER" -h -o "%i %j" | awk -v n="$NAME" '$2 ~ "^rb-" n "-(water|complex)$" {print $1}' | paste -sd:)
  if [ -n "$jobs" ]; then
    J=$(sbatch --parsable --dependency=afterany:$jobs -J "report-$NAME" "$HERE/report.sbatch" "$NAME")
    echo "not finished: ${missing[*]} -> resubmitted as $J after jobs $jobs"; exit 0
  fi
  fail "not finished and no queued job for: ${missing[*]}"
fi

# ---- 1. inputs, one torsion numbering ------------------------------------------------------------
quiet python -I "$HERE/prepare.py" --name "$NAME" || fail prepare
M=$WORK/manifest.json
read -r STRIDE MAXM <<< "$(python -I -c "import json; m=json.load(open('$M')); print(m['stride'], m['max_matches'] or '-')")"

# ---- 2. clustering, water and every complex ------------------------------------------------------
cluster() {   # cluster <system key in the manifest: water | complexes.<k>>
  local KEY=$1
  mapfile -t A < <(python -I - "$M" "$KEY" "$NAME" <<'EOF'
import json, os, sys
m = json.load(open(sys.argv[1])); key = sys.argv[2]
s = m["water"] if key == "water" else m["complexes"][int(key.split(".")[1])]
label = "water" if key == "water" else s["label"] + " complex"
print(s["label"])
args = ["--traj", s["traj"], "--top", s["top"], "--sdf", s["sdf"], "--stride", str(m["stride"]),
        "--no-single-cluster", "--name", f"{sys.argv[3]} {label}"]
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
  local LABEL=${A[0]}; local OUT=$WORK/$LABEL-clustering
  python "$SHUYU/cluster_prune.py" "${A[@]:1}" --rmsd-workers "${SLURM_CPUS_PER_TASK:-8}" --out "$OUT" > "$OUT.log" 2>&1
  local rc=$?
  [ $rc != 0 ] && [ $rc != 3 ] && { tail -20 "$OUT.log"; fail "cluster_prune $LABEL ($rc)"; }
  echo "   $LABEL: $(grep -hE '^n_clusters' "$OUT.log" | head -1) $( [ -s "$OUT/NO_CLUSTER.txt" ] && echo '(NO cluster found)')"
}
echo "== 2. clustering (no single cluster)"
cluster water
NC=$(python -I -c "import json; print(len(json.load(open('$M'))['complexes']))")
for ((k = 0; k < NC; k++)); do cluster "complexes.$k"; done

# ---- 3-5. complex comparisons --------------------------------------------------------------------
if [ "$NC" -gt 0 ]; then
  echo "== 3. RMSD density";       quiet python -I "$HERE/rmsd_density.py" --manifest "$M" --out "$WORK/rmsd_density" || fail rmsd_density
  echo "== 4. binding check";      python -I "$REPO/analysis/binding-check/check_bound.py" --thesis "$T" --runs "${CRUNS[@]}" \
                                      --name "$NAME" --title "$NAME in its protein: does the ligand stay in the pocket?" --out "$WORK" || fail check_bound
  echo "== 5. water vs complex";   quiet python -I "$HERE/water_vs_complex.py" --manifest "$M" --out "$WORK/water_vs_complex" || fail water_vs_complex
fi

# ---- results -------------------------------------------------------------------------------------
OUTDIR=$REPO/results/$NAME/$VERSION
mkdir -p "$OUTDIR"
for d in "$WORK"/*-clustering; do
  L=$(basename "$d" -clustering)
  for f in clusters_overview.png groups.csv NO_CLUSTER.txt; do [ -s "$d/$f" ] && cp "$d/$f" "$OUTDIR/${NAME}_${L}_${f/clusters_overview/clusters}"; done
done
if [ "$NC" -gt 0 ]; then
  cp "$WORK/rmsd_density.png" "$OUTDIR/${NAME}_rmsd_density.png"; cp "$WORK/rmsd_density.csv" "$OUTDIR/${NAME}_rmsd_density.csv"
  cp "$WORK/binding_$NAME.png" "$OUTDIR/${NAME}_binding.png";    cp "$WORK/binding_$NAME.csv" "$OUTDIR/${NAME}_binding.csv"
  cp "$WORK/water_vs_complex.png" "$OUTDIR/${NAME}_water_vs_complex.png"; cp "$WORK/water_vs_complex.csv" "$OUTDIR/${NAME}_water_vs_complex.csv"
fi
python -I - "$M" "$OUTDIR" "${SLURM_JOB_ID:-}" <<'EOF'
import json, os, sys, time, importlib.metadata as md_meta
m = json.load(open(sys.argv[1])); out = sys.argv[2]; T = "/scicore/home/lill0000/popa0002/thesis/"
tors = [l for k in ("torsions", "torsions_rigid", "torsions_polarH")
        for l in open(f"{m['water']['torsions']}/{k}.txt") if l.strip() and not l.startswith("#")]
L = [f"# {m['name']} -- {m['version']}", "",
     f"Made by `analysis/report/report.sbatch` (Slurm job {sys.argv[3]}, {time.strftime('%Y-%m-%d')}, "
     f"md-tools {md_meta.version('md-tools')}); state 0, stride {m['stride']}.", "",
     "| system | run | crystal reference |", "|---|---|---|",
     f"| water | `{os.path.relpath(m['water']['run'], T)}` ({' + '.join(os.path.basename(os.path.dirname(p)) for p in m['water']['parts'])}) | "
     + ", ".join(f"`{os.path.basename(r)}`" for r in m['water']['refs']) + " |"]
for c in m["complexes"]:
    L.append(f"| {c['label']} complex | `{c['run_rel']}` | `{os.path.basename(c['refs'][0])}` |")
L += ["", f"Torsions t0..t{len(tors) - 1} (the same in every figure; rotatable, then rigid, then polar-H; "
      "atom indices of the water ligand):", ""] + [f"- t{k}: `{l.split('#')[0].strip()}`  {('#' + l.split('#', 1)[1].strip()) if '#' in l else ''}" for k, l in enumerate(tors)]
L += ["", "Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), "
      "`*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: "
      "pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions)."]
open(f"{out}/README.md", "w").write("\n".join(L) + "\n")
EOF
echo "== $NAME done -> $OUTDIR  $(date)"
ls "$OUTDIR"
