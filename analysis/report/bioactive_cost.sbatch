#!/bin/bash
#SBATCH --job-name=bioactive-cost
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%j.out
# Cost of the bioactive conformation in water (dG_dist, dG_strain; Zivanovic et al. 2020) from the
# barrier-cut water clusters of report_cut.sbatch -> results/overview/bioactive-cost/. Extra arguments
# are passed on (e.g. --compounds folic-acid 4ibi).
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/bioactive_cost.py "$@"
