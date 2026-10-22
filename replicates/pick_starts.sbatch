#!/usr/bin/env bash
#SBATCH --job-name=pick-starts
#SBATCH --partition=scicore
#SBATCH --qos=1day
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/replicates/logs/%x-%A_%a.out
# Starting conformers for REST2 replicates (replicates/pick_starts.py), env analysis-env-main.
# Replicate 0 = the existing v4 water run, clustered with the same settings as run_symfirst.sbatch.
#   sbatch --array=0-1 replicates/pick_starts.sbatch      0 paracetamol v4   1 folate water v4
# Outputs: replicates/<name>/{candidates.csv,selection.json,starts/rep*.sdf,*.png}
set -euo pipefail
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$(/scicore/home/lill0000/popa0002/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
HERE=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/replicates
case ${SLURM_ARRAY_TASK_ID} in
0) NAME=paracetamol-v4; R=/scicore/home/lill0000/popa0002/thesis/PARA-REST2-v4
   ARGS=(--package-sdf "$R/build/parameter/molecule.sdf" --start-sdf "$R/build/parameter/TYL.sdf"
         --traj "$R/REST2-run1/solute_state0_prod1.nc" --top "$R/build/built.solute.pdb" --stride 1 --prune-rms 0.1) ;;
1) NAME=folic-acid-v4-water; R=/scicore/home/lill0000/popa0002/thesis/FOL-REST2-v4; T=$R/water/REST2-run1
   ARGS=(--package-sdf "$R/packages/FOL/molecule.sdf" --start-sdf "$R/water/build/FOL.sdf"
         --traj "$T/solute_state0_prod1.nc" --top "$R/water/build/built.solute.pdb" --stride 5
         --torsions-file "$T/torsions.txt" --torsions-file "$T/torsions_rigid.txt"
         --torsions-file "$T/torsions_polarH.txt" --no-single-cluster) ;;
esac
echo "== $NAME on $(hostname)  $(date)"
python "$HERE/pick_starts.py" "${ARGS[@]}" --n-pick 3 --name "$NAME" --out "$HERE/$NAME"
echo "== done  $(date)"
