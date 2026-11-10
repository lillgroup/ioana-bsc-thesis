#!/usr/bin/env bash
# Joint-clustering KL / JS: one array task per compound, then the figures once every task has finished.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
C=$(sbatch --parsable "$HERE/kl_joint.sbatch")
P=$(sbatch --parsable -J kl-joint-plot --array=0 --export=ALL,PLOT=1 --dependency=afterany:$C "$HERE/kl_joint.sbatch")
echo "joint $C | plot $P (after it)"
