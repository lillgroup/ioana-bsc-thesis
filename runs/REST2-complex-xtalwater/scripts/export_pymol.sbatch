#!/bin/bash
#SBATCH --job-name=xw-pymol
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/export-pymol-%j.out
# PyMOL movies (bsc-thesis/results/binding-check/export_pymol.py) for the 3QL3 runs: N1-H and 3H.
#   sbatch scripts/export_pymol.sbatch [system ...]     (default: FOL-3QL3 FOL-3QL3-3H)
R=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
for s in ${@:-FOL-3QL3 FOL-3QL3-3H}; do
  echo "######## $s"
  python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/results/binding-check/export_pymol.py \
    --run $R/$s --ligand FOL --out $R/$s/analysis/pymol
  ls -la $R/$s/analysis/pymol
done
