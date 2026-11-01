#!/usr/bin/env bash
# Submit REST2-rotbonds as a dependency chain (only sbatch calls; safe on the login node):
#   build (array task i)  ->  300-exchange tests, water and complex  ->  productions (60 / 40 ns)
# A test must succeed (afterok) before its production starts; a failed build or test leaves the rest
# pending with DependencyNeverSatisfied (scancel them, fix, rerun this for that compound).
#   bash scripts/submit.sh                 all 15
#   bash scripts/submit.sh 6cjr 4erf       some (names from compounds.csv)
#   NO_BUILD=1 bash scripts/submit.sh 6cjr # already built: tests + productions only
# GPU type: GPU=h200 (default a100) for partition, qos and gres. Job ids go to logs/jobs.txt.
set -euo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds
cd "$ROOT"
GPU=${GPU:-a100}
GPUARGS=(-p "$GPU" --ntasks=6 --gres="gpu:$GPU:3")
mapfile -t NAMES < <(tail -n +2 compounds.csv | cut -d, -f1)
want=("$@"); [ ${#want[@]} -gt 0 ] || want=("${NAMES[@]}")
idx=()
for w in "${want[@]}"; do
  k=-1; for j in "${!NAMES[@]}"; do [ "${NAMES[$j]}" = "$w" ] && k=$j; done
  [ $k -ge 0 ] || { echo "unknown compound $w"; exit 1; }
  idx+=("$k")
done
mkdir -p logs
if [ -z "${NO_BUILD:-}" ]; then
  BUILD=$(sbatch --parsable --array="$(IFS=,; echo "${idx[*]}")" scripts/build.sbatch)
  echo "build $BUILD (array ${idx[*]})" | tee -a logs/jobs.txt
fi
for k in "${idx[@]}"; do
  n=${NAMES[$k]}
  dep=(); [ -n "${BUILD:-}" ] && dep=(--dependency="afterok:${BUILD}_$k")
  for kind in water complex; do
    T=$(sbatch --parsable "${dep[@]}" -J "rb-test-$n-$kind" "${GPUARGS[@]}" --qos="$GPU-6hours" --time=06:00:00 \
        scripts/rest2.sbatch "$n-$kind/test/REST2-test")
    P=$(sbatch --parsable --dependency="afterok:$T" -J "rb-$n-$kind" "${GPUARGS[@]}" --qos="$GPU-1day" \
        scripts/rest2.sbatch "$n-$kind/REST2-run1")
    echo "$n-$kind test $T -> production $P ($GPU)" | tee -a logs/jobs.txt
  done
done
