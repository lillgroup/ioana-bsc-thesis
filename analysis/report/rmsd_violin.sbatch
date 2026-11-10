#!/bin/bash
#SBATCH --job-name=rmsd-violin
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%j.out
# RMSD-to-crystal violins, water vs complex (rmsd_violin.py), for the compounds given as arguments.
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
for n in "$@"; do python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/rmsd_violin.py --name "$n" 2>&1 | grep -iv warning; done
