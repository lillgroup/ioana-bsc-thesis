---
name: scicore-setup
description: The thesis stack was set up on sciCORE on 2026-10-22; paths, Slurm partitions/QOS, network limits and rebuild traps
metadata:
  type: project
---

On 2026-10-22 the stack was set up on sciCORE (Uni Basel), replacing sbrinz. Details are in SETUP.md §0.

- Envs `openmm-env-protonation` (b3603e0) and `analysis-env-main` (69c262b) live in
  `/scicore/home/lill0000/popa0002/software/md-stack`. The simulation env matches the sbrinz lock.
  The analysis env was built from environment.yml, not cloned, so it is a superset with minor drift
  (rdkit 2026.03.1, numpy 2.5.4, mdtraj 1.11.2).
- `MD_DATA=/scicore/home/lill0000/popa0002/md-data`. The paracetamol package is in
  `parameters/ligands/CHEMBL112/param_e932f4c4f371`. `user.config` keeps `person_id: ipopa`.
- All sbatch headers were ported: GPU jobs use a100 + a100-1day, CPU jobs use scicore + 1day.
  The 5-GPU folate layouts became 4 GPUs (3/3/2/2).
- The GPU test passed: job 28226867, CUDA, 906 ns/day per replica.

**Why:** the move off sbrinz, done after the user asked for "prepare MD REST2 so we can run new simulations".
**How to apply:**
- Before any rebuild, move `build/ligands/` aside and use `build-top --overwrite`.
- Never run the `--parameterize` line of the v4 build scripts. See [[slurm-always]] and [[md-tools-openmm-env]].
