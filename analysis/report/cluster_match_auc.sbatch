#!/bin/bash
#SBATCH --job-name=cluster-match
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-19
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%A_%a.out
# Matched-water-population curve and its area score (cluster_match_auc.py), one task per compound with a
# complex; PLOT=1 draws the overview from the per-compound CSVs. --array=0 = folic acid.
set -uo pipefail
mapfile -t COMPOUNDS < <(awk -F, 'NR > 1 && $5 != "" {print $1}' /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/compounds.csv)
K=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/cluster_match_auc.py
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=1
if [ "${PLOT:-0}" = 1 ]; then micromamba activate analysis-env-main; exec python -I $K --plot; fi
micromamba activate analysis-env-cut
python -I $K --name "${COMPOUNDS[$SLURM_ARRAY_TASK_ID]}" --workers "${SLURM_CPUS_PER_TASK:-8}" 2>&1 | grep -iv warning
exit "${PIPESTATUS[0]}"
