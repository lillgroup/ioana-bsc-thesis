#!/bin/bash
#SBATCH --job-name=xw-frag
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/fragments-%j.out
# Folate per-fragment drift (pteridine / pABA / glutamate) and 3QL3 contacts, with the v4 script.
R=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
F=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/runs/FOL-REST2-v4/scripts/fragments.py
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
python -I $F $R/FOL-3QL3/REST2-run1 $R/FOL-3QL3/build/built.solute.pdb contacts 2>&1 | grep -iv warning | tee $R/FOL-3QL3/analysis/fragments.txt
python -I $F $R/FOL-4LRH/REST2-run1 $R/FOL-4LRH/build/built.solute.pdb 2>&1 | grep -iv warning | tee $R/FOL-4LRH/analysis/fragments.txt
