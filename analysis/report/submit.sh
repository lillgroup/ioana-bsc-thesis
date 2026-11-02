#!/usr/bin/env bash
# One report.sbatch job per compound of compounds.csv (only sbatch calls; safe on the login node).
# A compound whose runs are still going waits by itself (report.sbatch resubmits after its run jobs).
#   bash analysis/report/submit.sh               every compound
#   bash analysis/report/submit.sh 4ibi 6cjr     some
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
mapfile -t ALL < <(tail -n +2 "$HERE/compounds.csv" | cut -d, -f1)
for n in "${@:-${ALL[@]}}"; do
  printf '%s ' "${ALL[@]}" | grep -qw -- "$n" || { echo "unknown compound $n"; exit 1; }
  echo "$n $(sbatch --parsable -J "report-$n" "$HERE/report.sbatch" "$n")"
done
