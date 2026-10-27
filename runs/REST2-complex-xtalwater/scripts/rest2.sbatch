#!/bin/bash
#SBATCH --job-name=xw-rest2
#SBATCH --partition=a100
#SBATCH --qos=a100-1day
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:a100:4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/rest2-%x-%j.out
#SBATCH --signal=B:USR1@900
# One complex REST2 run with crystal waters (min, eq_1..3, production) inside its Slurm GPU allocation.
# Unattended: 15 min before the time limit, or when run.sh fails, the job resubmits itself with the
# same resources (try N of 4; the optional 2nd argument is the try number), and run.sh resumes from
# its checkpoint. A completed run writes check_run.py's acceptance/binding report into this log.
#   folate (10 replicas, 4 GPUs = 3/3/2/2 ranks):  sbatch -J xw-FOL-3QL3 scripts/rest2.sbatch FOL-3QL3/REST2-run1
#   small molecules (6 replicas, 3 GPUs):          sbatch -J xw-4ibi --ntasks=6 --gres=gpu:a100:3 scripts/rest2.sbatch 4ibi/REST2-run1
#
# 10 replicas on 4 GPUs (sciCORE a100 nodes have 4) = 3/3/2/2 ranks per GPU, which md-tools allows only under NVIDIA MPS. The MPS
# daemon is private to this job (pipe dir keyed by job id) and is started on exactly the GPUs Slurm
# allocated; clients then see them as 0..n-1. It is stopped when the job ends.
# run.sh resumes an interrupted stage from its checkpoint, so resubmitting after a crash is safe.
set -uo pipefail
ROOT=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
RUN="$ROOT/$1"
[ -x "$RUN/run.sh" ] || { echo "no $RUN/run.sh"; exit 2; }
TRY=${2:-1}; MAX_TRIES=4
# this job's own request, read back so a resubmission asks for exactly the same
NT=$SLURM_NTASKS; PART=$SLURM_JOB_PARTITION
QOS=$(squeue -h -j "$SLURM_JOB_ID" -o %q); TL=$(squeue -h -j "$SLURM_JOB_ID" -o %l)
GRES=$(squeue -h -j "$SLURM_JOB_ID" -o %b); GRES=${GRES#gres/}; GRES=${GRES#gres:}
resubmit() {
  [ -n "${RESUBMITTED:-}" ] && return; RESUBMITTED=1
  if [ "$TRY" -ge "$MAX_TRIES" ]; then echo "NOT resubmitting: try $TRY of $MAX_TRIES"; return; fi
  # a clean login environment, so none of this job's SLURM_* variables leak into the new one
  env -i HOME="$HOME" PATH=/usr/bin:/bin bash -lc "sbatch -J '$SLURM_JOB_NAME' -p '$PART' --qos='$QOS' \
    --ntasks=$NT --gres='$GRES' --time='$TL' '$ROOT/scripts/rest2.sbatch' '$1' $((TRY + 1))" \
    && echo "resubmitted as try $((TRY + 1)) ($1)  $(date)"
}
trap 'echo "time limit in 15 min  $(date)"; resubmit "$1"' USR1
echo "try $TRY of $MAX_TRIES  ($PART, $QOS, $NT tasks, $GRES, $TL)"

export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$(/scicore/home/lill0000/popa0002/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation   # Shu-Yu branch b3603e0 (ligand pH model, REST2 v4)
export MD_DATA=/scicore/home/lill0000/popa0002/md-data

ALLOC="${CUDA_VISIBLE_DEVICES:?Slurm set no CUDA_VISIBLE_DEVICES}"
NG=$(awk -F, '{print NF}' <<<"$ALLOC")
echo "host $(hostname)  job $SLURM_JOB_ID  GPUs (Slurm) $ALLOC  run $1  $(date)"
nvidia-smi -L

export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps-$SLURM_JOB_ID/pipe"
export CUDA_MPS_LOG_DIRECTORY="/tmp/mps-$SLURM_JOB_ID/log"
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
stop_mps() { echo quit | nvidia-cuda-mps-control 2>/dev/null; rm -rf "/tmp/mps-$SLURM_JOB_ID"; }
trap stop_mps EXIT
export CUDA_DEVICE_ORDER=PCI_BUS_ID
nvidia-cuda-mps-control -d || { echo "could not start MPS"; exit 1; }
echo "MPS daemon on $ALLOC ($CUDA_MPS_PIPE_DIRECTORY)"

export CUDA_VISIBLE_DEVICES=$(seq -s, 0 $((NG - 1)))   # the daemon's devices, as clients see them
# md-tools takes the MPI world size from OMPI_COMM_WORLD_SIZE, else PMI_SIZE, else SLURM_NTASKS. The
# serial stages (min, eq_1..3) are not under mpirun, so SLURM_NTASKS=10 made them refuse ("launcher
# world size 10, MPI communicator size 1"; jobs 8024/8025). mpirun still takes its 10 slots from
# SLURM_TASKS_PER_NODE, and its ranks get OMPI_COMM_WORLD_SIZE, which is read first.
unset SLURM_NTASKS SLURM_NPROCS
cd "$RUN"
# in the background, so the USR1 trap runs while run.sh is still going
./run.sh & RS=$!
while :; do wait $RS; rc=$?; kill -0 $RS 2>/dev/null || break; done
echo "run.sh rc=$rc  $(date)"
grep -E 'throughput|run_status' remd_records/REST2_prod1.out 2>/dev/null
if grep -q 'run_status: completed' remd_records/REST2_prod1.out 2>/dev/null; then
  LIG=$(awk '/ligand_scaling_dict:/{getline; gsub(/[ :]/, ""); print; exit}' ../build/scaler.config)
  TOP=$(( $(grep -oE 'number_of_replicas *= *[0-9]+' ../input/REST2.in 2>/dev/null | grep -oE '[0-9]+$' || echo 1) - 1 ))
  echo "== check ($LIG, states 0,$TOP)"
  python -I "$ROOT/scripts/check_run.py" "$RUN" --states "0,$TOP" --ligand "$LIG" 2>&1 | grep -v Warning
  exit 0
fi
resubmit "$1"
exit $rc
