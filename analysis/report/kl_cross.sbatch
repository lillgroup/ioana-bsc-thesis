#!/bin/bash
#SBATCH --job-name=kl-cross
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --array=0-19
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/analysis-logs/%x-%A_%a.out
# KL in both directions with Shu-Yu's protocol (clusters of one ensemble, frames of the other assigned),
# on the existing state-0 ensembles -- no MD. One array task per compound (COMPOUNDS below), method from
# METHOD (cut: env analysis-env-cut, 306b109; hdbscan: env analysis-env-main, c30c1c0). PLOT=1 instead
# draws both figures from the per-compound CSVs. bash kl_cross_submit.sh does all of it.
#   sbatch --export=ALL,METHOD=cut kl_cross.sbatch          sbatch --array=0 --export=ALL,PLOT=1 kl_cross.sbatch
set -uo pipefail
# every compound of compounds.csv with a complex (paracetamol has none), in file order
mapfile -t COMPOUNDS < <(awk -F, 'NR > 1 && $5 != "" {print $1}' /scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/compounds.csv)
K=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/analysis/report/kl_cross.py
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
if [ "${PLOT:-0}" = 1 ]; then micromamba activate analysis-env-main; exec python -I $K --plot; fi
case ${METHOD:?METHOD=cut or hdbscan} in cut) micromamba activate analysis-env-cut ;; hdbscan) micromamba activate analysis-env-main ;; esac
python -I $K --method "$METHOD" --compounds "${COMPOUNDS[$SLURM_ARRAY_TASK_ID]}" 2>&1 | grep -iv warning
exit "${PIPESTATUS[0]}"
