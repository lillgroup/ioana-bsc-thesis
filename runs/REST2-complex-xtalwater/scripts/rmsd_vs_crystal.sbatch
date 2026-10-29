#!/bin/bash
#SBATCH --job-name=xw-rmsd
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/rmsd-vs-crystal-%j.out
# Folate RMSD to the crystal conformation (folate superposed on folate) for FOL-4LRH and FOL-3QL3 with
# crystal waters, and the density plot against the earlier v4 runs without them (read from the CSV
# in results/symfirst-2026-10-15/1_v4-new-rest2/). See scripts/rmsd_vs_crystal.py.
#   sbatch scripts/rmsd_vs_crystal.sbatch      (from REST2-complex-xtalwater/)
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
T=/scicore/home/lill0000/popa0002/thesis
mkdir -p $T/REST2-complex-xtalwater/analysis/rmsd-vs-crystal
python -I $T/bsc-thesis/runs/REST2-complex-xtalwater/scripts/rmsd_vs_crystal.py --thesis $T \
  --out $T/REST2-complex-xtalwater/analysis/rmsd-vs-crystal/folic-acid-v4-xtalwater_rmsd_vs_crystal \
  2>&1 | grep -iv warning
