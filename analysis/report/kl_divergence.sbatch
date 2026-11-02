#!/bin/bash
#SBATCH --job-name=kl-divergence
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%j.out
# KL(complex || water) on one shared clustering per (ligand, complex), both methods, each in its own
# environment (barrier cuts: analysis-env-cut, 306b109; HDBSCAN: analysis-env-main, c30c1c0), then the
# two figures -> results/overview/kl-divergence/. Extra arguments go to both runs (e.g. --compounds 4ibi).
set -uo pipefail
K=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/kl_divergence.py
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
micromamba activate analysis-env-cut;  python -I $K --method cut "$@" 2>&1 | grep -iv warning || exit 1
micromamba activate analysis-env-main; python -I $K --method hdbscan "$@" 2>&1 | grep -iv warning || exit 1
python -I $K --plot
