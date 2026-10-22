#!/bin/bash
#SBATCH --job-name=env-verify
#SBATCH --partition=a100
#SBATCH --qos=a100-30min
#SBATCH --gres=gpu:a100:1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/scicore/home/lill0000/popa0002/software/logs/%x-%j.out
# SETUP.md §1.3 checks on a GPU node, plus the §1.4 machine config (written once; skipped if present).
set -uo pipefail
SOFT=/scicore/home/lill0000/popa0002/software
export MAMBA_ROOT_PREFIX=$SOFT/md-stack
eval "$($SOFT/md-stack/bin/micromamba shell hook -s bash)"
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
nvidia-smi -L
micromamba activate openmm-env-protonation
md-openmm --version; command -v sqm mpiexec nvidia-cuda-mps-control
python -c "from openmm import Platform as P; print([P.getPlatform(i).getName() for i in range(P.getNumPlatforms())])"
python -c "
import openmm as mm, openmm.unit as u
s = mm.System(); s.addParticle(1.0)
c = mm.Context(s, mm.VerletIntegrator(1*u.femtosecond), mm.Platform.getPlatformByName('CUDA'), {'Precision': 'mixed'})
print('CUDA context OK on', c.getPlatform().getPropertyValue(c, 'DeviceName'))"
python -c "from md_tools.openmm.ligand_protonation import openbabel_version as v; print('openbabel', v())"
CFG=${XDG_CONFIG_HOME:-$HOME/.config}/md-tools/user.config
[ -f "$CFG" ] || md-openmm data-register --init --noninteractive --name ipopa --person-id ipopa --md-data "$MD_DATA"
cat "$CFG"
micromamba activate analysis-env-main
python -c "from md_tools.analysis import t_hdbscan, THDBSCAN; print('align_all', hasattr(THDBSCAN,'align_all'))"
python -c "import sklearn, networkx, matplotlib; print('sklearn', sklearn.__version__, 'networkx', networkx.__version__)"
