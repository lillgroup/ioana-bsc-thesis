#!/bin/bash
#SBATCH --job-name=entropy-scan
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%j.out
# S_c vs barrier-cut settings (pi_min x E*), entropy_scan.py, for the compounds given as arguments.
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-cut
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
for n in "$@"; do python -I /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/entropy_scan.py --name "$n" 2>&1 | grep -iv warning; done
