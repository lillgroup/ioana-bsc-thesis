#!/bin/bash
#SBATCH --job-name=xw-check
#SBATCH --partition=scicore
#SBATCH --qos=30min
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater/logs/check-%j.out
# Acceptance + binding of every finished run: <system>/REST2-run1, and <system>/test/REST2-test.
R=/scicore/home/lill0000/popa0002/thesis/REST2-complex-xtalwater
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
for x in FOL-3QL3:FOL:9 FOL-4LRH:FOL:9 4ibi:1D8:5 4ibk:1DE:5 4j7d:I31:5 5mwh:UWX:5 5ov8:AXN:5; do
  IFS=: read -r s lig top <<<"$x"
  for run in REST2-run1 test/REST2-test; do
    grep -q 'run_status: completed' "$R/$s/$run/remd_records/REST2_prod1.out" 2>/dev/null || continue
    echo "######## $s $run"
    grep -hE 'throughput' "$R/$s/$run/remd_records/REST2_prod1.out"
    python -I "$R/scripts/check_run.py" "$R/$s/$run" --states 0,$top --ligand "$lig" 2>&1 | grep -v Warning
  done
done
