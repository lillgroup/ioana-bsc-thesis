#!/usr/bin/env bash
#SBATCH --job-name=ana-folpara
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=0-1
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/results/slurm-%x-%A_%a.out
# Latest Shu-Yu analysis (cluster_prune.py + cluster_convergence.py) for
#   0: folic acid, FOLn-REST2-50ns (100 ns, 8 replicas), crystal refs 4LRH chain A + 3QL3
#   1: paracetamol, REST2-PARA (run interrupted at 3.5 ns; stride 1, 0.5 ns windows, no crystal ref)
set -euo pipefail

export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$(/scicore/home/lill0000/popa0002/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

SHUYU=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/shuyu
CRYSTAL=/scicore/home/lill0000/popa0002/thesis/bsc-thesis/crystal

case ${SLURM_ARRAY_TASK_ID} in
0)  NAME=folic-acid; ROOT=/scicore/home/lill0000/popa0002/thesis/FOLn-REST2-50ns
    ARGS=(--sdf build/FOL.sdf
          --torsions-file REST2-run1/torsions.txt --torsions-file REST2-run1/torsions_rigid.txt
          --torsions-file REST2-run1/torsions_polarH.txt --stride 5)
    REFS=(--ref-sdf "$CRYSTAL/fol_4LRH_A.sdf" --ref-sdf "$CRYSTAL/fol_3QL3_A.sdf")
    STEP=5 ;;
1)  NAME=paracetamol; ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-PARA
    # no --torsions-file: cluster_prune.py's default set is the paracetamol tutorial's
    # omega / aryl / phenol torsions
    ARGS=(--sdf build/TYL.sdf --stride 1)
    REFS=()
    STEP=0.5 ;;
esac

cd "$ROOT"
COMMON=(--traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb "${ARGS[@]}"
        --no-single-cluster --name "$NAME" --rest2-config REST2.config)

echo "== $NAME cluster_prune.py  $(date)"
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" \
    --out "$NAME-clustering-shuyu-nosingle-ref" > "$NAME-clustering-shuyu-nosingle-ref.log" 2>&1

echo "== $NAME cluster_convergence.py  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns "$STEP" --plateau 3 \
    --out "$NAME-convergence" > "$NAME-convergence.log" 2>&1

echo "== $NAME done  $(date)"
grep -h "converged" "$NAME-convergence.log" | tail -1
