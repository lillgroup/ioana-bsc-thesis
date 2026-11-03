# 4crf -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365393, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4crf-water` (REST2-run1) | `4crf.input.sdf` |
| 4CRF complex | `REST2-rotbonds/4crf-complex` | `4crf.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `3 2 0 1`  #X141-C3 X141-C2 X141-C1 X141-O1 | rotatable
- t1: `16 15 14 0`  #X141-C12 X141-C11 X141-N2 X141-C1 | rotatable
- t2: `17 16 15 14`  #X141-O4 X141-C12 X141-C11 X141-N2 | rotatable
- t3: `19 18 15 14`  #X141-C14 X141-C13 X141-C11 X141-N2 | rotatable
- t4: `20 19 18 15`  #X141-C15 X141-C14 X141-C13 X141-C11 | rotatable
- t5: `38 37 35 30`  #X141-N4 X141-C29 X141-C28 X141-C24 | rotatable
- t6: `16 38 37 35`  #X141-C12 X141-N4 X141-C29 X141-C28 | rotatable
- t7: `15 14 0 1`  #X141-C11 X141-N2 X141-C1 X141-O1 | rigid amide-like C-N
- t8: `37 38 16 15`  #X141-C29 X141-N4 X141-C12 X141-C11 | rigid amide-like C-N
- t9: `6 8 13 44`  #X141-C6 X141-C8 X141-O3 X141-H6 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29468045, 2026-11-03);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 150 clusters | 13.44, 6.17, 5.72, 5.52, 5.49, 3.65, 3.59, 2.92, 2.84, 2.79 ... (+140 more, 47.8 % together) |
| 4CRF | 2 clusters | 89.95, 10.05 |
