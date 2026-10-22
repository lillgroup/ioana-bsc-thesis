#!/bin/bash
#SBATCH --job-name=env-protonation
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scicore/home/lill0000/popa0002/software/logs/%x-%j.out
# Simulation env (SETUP.md §1.2): bootstraps micromamba + the MD-tools clone and its two pinned
# worktrees if missing, then env openmm-env-protonation, editable from MD-tools-protonation (b3603e0).
# sciCORE: compute nodes cannot reach github.com or the micromamba download (403 via proxy), so run
# the curl and the git clone/worktree lines once on the login node first; conda-forge works from jobs.
set -euo pipefail
SOFT=/scicore/home/lill0000/popa0002/software
mkdir -p $SOFT/md-stack
[ -x $SOFT/md-stack/bin/micromamba ] || \
  (cd $SOFT/md-stack && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba)
export MAMBA_ROOT_PREFIX=$SOFT/md-stack
eval "$($SOFT/md-stack/bin/micromamba shell hook -s bash)"
[ -d $SOFT/MD-tools/.git ] || git clone https://github.com/csy0000/MD-tools.git $SOFT/MD-tools
[ -d $SOFT/MD-tools-protonation ] || git -C $SOFT/MD-tools worktree add $SOFT/MD-tools-protonation b3603e0
[ -d $SOFT/MD-tools-main ]        || git -C $SOFT/MD-tools worktree add $SOFT/MD-tools-main 69c262b
SRC=$SOFT/MD-tools-protonation
echo "== $(git -C $SRC log --oneline -1)  on $(hostname)  $(date)"
micromamba create -y -n openmm-env-protonation -f "$SRC/environment.yml"
micromamba activate openmm-env-protonation
pip install --no-deps -e "$SRC"
echo "== check"
python -c "import md_tools, importlib.metadata as m, os; print('md-tools', m.version('md-tools'), 'from', os.path.dirname(md_tools.__file__))"
python -c "from md_tools.openmm.ligand_protonation import openbabel_version; print('openbabel', openbabel_version())"
python -c "import openmm; print('openmm', openmm.__version__)"
which md-openmm
echo "== done $(date)"
