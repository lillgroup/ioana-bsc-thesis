#!/bin/bash
#SBATCH --job-name=conf-entropy
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%j.out
# Conformational entropy e^S_c from the cluster populations of report.sbatch (HDBSCAN) and
# report_cut.sbatch (barrier cuts) -> results/overview/conformational-entropy/. Extra arguments are
# passed on (e.g. --compounds folic-acid 4ibi).
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/conformational_entropy.py "$@"
