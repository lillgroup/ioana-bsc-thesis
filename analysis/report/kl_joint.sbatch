#!/bin/bash
#SBATCH --job-name=kl-joint
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-19
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%A_%a.out
# KL and Jensen-Shannon divergence on ONE joint barrier-cut clustering of the water and complex frames
# (kl_joint.py), on the existing state-0 ensembles -- no MD. One array task per compound with a complex
# (env analysis-env-cut); PLOT=1 instead draws the figures (env analysis-env-main). bash kl_joint_submit.sh
# does both.
set -uo pipefail
mapfile -t COMPOUNDS < <(awk -F, 'NR > 1 && $5 != "" {print $1}' /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/compounds.csv)
K=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/kl_joint.py
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
if [ "${PLOT:-0}" = 1 ]; then micromamba activate analysis-env-main; exec python -I $K --plot; fi
micromamba activate analysis-env-cut
python -I $K --compounds "${COMPOUNDS[$SLURM_ARRAY_TASK_ID]}" 2>&1 | grep -iv warning
exit "${PIPESTATUS[0]}"
