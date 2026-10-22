#!/bin/bash
#SBATCH --job-name=gpu-test-para
#SBATCH --partition=a100
#SBATCH --qos=a100-6hours
#SBATCH --gres=gpu:a100:2
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scicore/home/lill0000/popa0002/software/logs/%x-%j.out
# SETUP.md §6 "one short GPU test" on a new machine: PARA-REST2-v4 (4 replicas, 2 ranks per GPU under
# MPS) copied to a throwaway folder, built from the COMMITTED package (no --parameterize: the line in
# runs/PARA-REST2-v4/scripts/build.sbatch is skipped on purpose), number_of_exchanges cut to 300.
# Passes when built.out says "reused", the records say CUDA, and every neighbour pair exchanges.
set -euo pipefail
THESIS=/scicore/home/lill0000/popa0002/thesis
TEST=$THESIS/tmp-gpu-test/PARA-REST2-v4-$SLURM_JOB_ID
SOFT=/scicore/home/lill0000/popa0002/software
[ -e "$TEST" ] && { echo "$TEST exists; remove it or pick another folder"; exit 2; }
mkdir -p "$(dirname "$TEST")"
cp -a "$THESIS/bsc-thesis/runs/PARA-REST2-v4" "$TEST"
sed -i 's/^\(  number_of_exchanges:\) *[0-9]*/\1 300/' "$TEST/REST2.config"
grep number_of_exchanges "$TEST/REST2.config"

export MAMBA_ROOT_PREFIX=$SOFT/md-stack
eval "$($SOFT/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
echo "== $(md-openmm --version)  on $(hostname)  $(date)"

cd "$TEST/build"
cp TYL.sdf TYL.committed.sdf   # committed output of the sbrinz build; rebuilt below and compared
# build/ligands/ in git is a partial build output (no parameters.ffxml); build-top refuses it, so move it aside
mv ligands ligands.committed
md-openmm build-top --overwrite -i parameter/TYL.sdf --config build-top.config -os built.xml -op built.pdb -log built.log > built.out 2>&1
cmp -s TYL.sdf TYL.committed.sdf && echo "TYL.sdf identical to committed" || { diff TYL.committed.sdf TYL.sdf | head -20 || true; }
grep -E 'reused|ligand protonation|solvation|atoms  |status' built.out || true
md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1
grep -E 'selection|unscaled torsions|ligand amides|state[0-9]+\.xml' scaler.out || true
cd "$TEST"
md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1
grep -E 'status' build-md.out || true; grep -hE 'number_of_replicas|number_of_exchanges' input/REST2.in || true

ALLOC="${CUDA_VISIBLE_DEVICES:?Slurm set no CUDA_VISIBLE_DEVICES}"
NG=$(awk -F, '{print NF}' <<<"$ALLOC")
nvidia-smi -L
export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps-$SLURM_JOB_ID/pipe"
export CUDA_MPS_LOG_DIRECTORY="/tmp/mps-$SLURM_JOB_ID/log"
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
stop_mps() { echo quit | nvidia-cuda-mps-control 2>/dev/null; rm -rf "/tmp/mps-$SLURM_JOB_ID"; }
trap stop_mps EXIT
export CUDA_DEVICE_ORDER=PCI_BUS_ID
nvidia-cuda-mps-control -d || { echo "could not start MPS"; exit 1; }
export CUDA_VISIBLE_DEVICES=$(seq -s, 0 $((NG - 1)))
unset SLURM_NTASKS SLURM_NPROCS
cd "$TEST/REST2-run1"
set +e
./run.sh
rc=$?
set -e
echo "run.sh rc=$rc  $(date)"
grep -E 'throughput|run_status' remd_records/REST2_prod1.out 2>/dev/null
grep -rhoE '"resolved_platform": *"[A-Za-z]+"' remd_records | sort | uniq -c
python -I - <<'PY'
import json
cr = json.load(open("remd_records/restart_prod1.json")).get("completion_report", {})
for p in cr.get("by_neighbouring_pair", []):
    print("pair", p)
PY
exit $rc
