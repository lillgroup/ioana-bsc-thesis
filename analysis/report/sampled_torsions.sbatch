#!/bin/bash
#SBATCH --job-name=sampled-torsions
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%j.out
# sampled_torsions.py for one system:  sbatch analysis/report/sampled_torsions.sbatch folic-acid 3QL3
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-cut
python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/sampled_torsions.py --name "$1" --system "$2" 2>&1 | grep -iv warning
