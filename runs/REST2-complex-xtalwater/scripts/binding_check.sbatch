#!/bin/bash
#SBATCH --job-name=xw-binding
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/binding-check-%j.out
# bsc-thesis/results/binding-check/check_bound.py on the crystal-water runs (the older runs are not on
# this cluster and are skipped). Pocket-frame RMSD, ligand-pocket COM, contacts, native Q; state 0.
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/results/binding-check/check_bound.py \
  --thesis /scicore/home/lill0000/popa0002/thesis \
  --out /scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/analysis/binding-check
