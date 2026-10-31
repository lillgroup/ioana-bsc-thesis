#!/bin/bash
#SBATCH --job-name=xw-ana-3H
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/analyse-3h-%j.out
# Tautomer test analysis: per-fragment drift + Asp27/Arg57/NADP+ contacts for FOL-3QL3-3H, and the
# binding check with all crystal-water runs incl. FOL-3QL3-3H (separate folder: binding-check-3H).
R=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
B=/scicore/home/lill0000/popa0002/thesis/bsc-thesis
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
echo "######## fragments FOL-3QL3-3H"
python -I $B/runs/FOL-REST2-v4/scripts/fragments.py $R/FOL-3QL3-3H/REST2-run1 $R/FOL-3QL3-3H/build/built.solute.pdb contacts \
  2>&1 | grep -iv warning | tee $R/FOL-3QL3-3H/analysis/fragments.txt
micromamba activate analysis-env-main
echo "######## binding check"
python -I $B/results/binding-check/check_bound.py --thesis /scicore/home/lill0000/popa0002/thesis --out $R/analysis/binding-check-3H 2>&1 | grep -v '^skip'
