# Shared env for the FOL REST2 workflow
export MAMBA_ROOT_PREFIX=/scicore/home/lill0000/popa0002/software/md-stack
eval "$(/scicore/home/lill0000/popa0002/software/md-stack/bin/micromamba shell hook -s bash)" 2>/dev/null
micromamba activate openmm-env
export SCRIPTS=/scicore/home/lill0000/popa0002/thesis/bsc-thesis
export ROOT=/scicore/home/lill0000/popa0002/thesis/FOLn-REST2
export RES=FOL
export MD_DATA=/scicore/home/lill0000/popa0002/md-data
