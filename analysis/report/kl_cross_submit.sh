#!/usr/bin/env bash
# Both methods (one array task per compound each), then the figures once every task has finished.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
C=$(sbatch --parsable -J kl-cross-cut --export=ALL,METHOD=cut "$HERE/kl_cross.sbatch")
H=$(sbatch --parsable -J kl-cross-hdbscan --export=ALL,METHOD=hdbscan --cpus-per-task=8 "$HERE/kl_cross.sbatch")
P=$(sbatch --parsable -J kl-cross-plot --array=0 --export=ALL,PLOT=1 --dependency=afterany:$C:$H "$HERE/kl_cross.sbatch")
echo "cut $C | hdbscan $H | plot $P (after both)"
