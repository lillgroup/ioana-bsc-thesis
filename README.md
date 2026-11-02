# bsc-thesis — REST2 of drug-like ligands in water and bound, with torsion clustering

REST2 simulations (MD-tools / OpenMM) of folate, paracetamol and small drug-like ligands with 4-7
rotatable bonds, each in water and in its protein, analysed with Shu-Yu's symmetry-first torsion
clustering (MD-tools `t_hdbscan`).

| folder | contents |
|---|---|
| `results/` | **results only**, one folder per compound, one subfolder per simulation version (see `results/README.md`); `results/overview/` for figures that cover several compounds |
| `analysis/` | **analysis code**: `shuyu_clustering/` (torsion selection, clustering, convergence, bioactive RMSD, water-vs-complex), `binding-check/` (pocket-frame RMSD, COM, contacts, PyMOL export), `run_*.sbatch` / `run_*.sh` drivers |
| `ligand-selection/` | choosing ligands and conformers by rotatable bonds: `filter.py`, `easy_candidates.csv`, `rotbonds_histogram.png`, `conformer_generation_from_SMILES.py`, `conformers/`; `rotbonds/` = the 4-7-bond PDBbind refined-set scan and structure screen |
| `runs/` | **simulation inputs only**, one folder per dataset: configs, parameter packages, prepared complexes, build / run scripts (see `runs/README.md`) |
| `crystal/` | crystal ligand SDFs used as references |
| `envs/` | exact environment specs and install / test jobs |
| `docs/notes/` | project notes (the why behind protocol choices) |

No simulation outputs in git (trajectories, checkpoints, built systems, run logs): they stay in the
run folders on the cluster (`$THESIS/<dataset>/`). Moving to a new machine: `SETUP.md`.
