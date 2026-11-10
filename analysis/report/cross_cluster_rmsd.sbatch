#!/bin/bash
#SBATCH --job-name=cross-rmsd
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-19
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%A_%a.out
# Cluster-by-cluster RMSD, complex clusters x water clusters (cross_cluster_rmsd.py), one array task per
# compound with a complex, on the existing barrier-cut runs. --array=0 = folic acid only.
set -uo pipefail
mapfile -t COMPOUNDS < <(awk -F, 'NR > 1 && $5 != "" {print $1}' /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/compounds.csv)
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-cut
export OMP_NUM_THREADS=1
python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/cross_cluster_rmsd.py \
  --name "${COMPOUNDS[$SLURM_ARRAY_TASK_ID]}" --workers "${SLURM_CPUS_PER_TASK:-8}" 2>&1 | grep -iv warning
exit "${PIPESTATUS[0]}"
