#!/bin/bash
#SBATCH --job-name=env-analysis-main
#SBATCH --partition=scicore
#SBATCH --qos=6hours
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scicore/home/lill0000/popa0002/software/logs/%x-%j.out
# Analysis env (SETUP.md §1.2): the MD-tools stack + md-tools main 69c262b (pinned, non-editable)
# + scikit-learn/networkx/matplotlib. Needs the clone/worktrees from install_protonation_env.sbatch.
set -euo pipefail
SOFT=/scicore/home/lill0000/popa0002/software
export MAMBA_ROOT_PREFIX=$SOFT/md-stack
eval "$($SOFT/md-stack/bin/micromamba shell hook -s bash)"
SRC=$SOFT/MD-tools-main
echo "== $(git -C $SRC log --oneline -1)  on $(hostname)  $(date)"
micromamba create -y -n analysis-env-main -f "$SRC/environment.yml"
micromamba activate analysis-env-main
micromamba install -y -c conda-forge scikit-learn networkx matplotlib
pip install --no-deps "$SRC"
python - <<'PY'
import md_tools, importlib.metadata as m, os
from md_tools.analysis import t_hdbscan, THDBSCAN, TorsionDefinitions
print("md-tools", m.version("md-tools"), "from", os.path.dirname(md_tools.__file__))
print("align_all:", hasattr(THDBSCAN, "align_all"), " density_peak:", hasattr(THDBSCAN, "density_peak"))
PY
echo "== done $(date)"
