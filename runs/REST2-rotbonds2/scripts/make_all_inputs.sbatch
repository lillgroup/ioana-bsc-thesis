#!/bin/bash
#SBATCH --job-name=rb2-inputs
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds2/logs/inputs-%j.out
# make_inputs.py, make_configs.py, prep.py for every compound of compounds.csv (or the names given),
# one at a time, so one failure does not stop the others.
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-rotbonds2
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-reduce
cd "$ROOT"
names=("$@"); [ ${#names[@]} -gt 0 ] || mapfile -t names < <(tail -n +2 compounds.csv | cut -d, -f1)
for n in "${names[@]}"; do
  echo "=== $n"
  python -I scripts/make_inputs.py "$ROOT" "$n" 2>&1 | tail -3 || true
  python -I scripts/make_configs.py "$ROOT" "$n" 2>&1 | tail -1
  python -I scripts/prep.py "$ROOT" "$n" 2>&1 | tail -6
done
