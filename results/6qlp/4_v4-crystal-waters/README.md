# 6qlp -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620728, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/6qlp-water` (REST2-run1) | `6qlp.input.sdf` |
| 6QLP complex | `REST2-rotbonds/6qlp-complex` | `6qlp.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `13 12 0 28`  #X171-C14 X171-C13 X171-C1 X171-O2 | rotatable
- t1: `5 16 1 15`  #X171-C6 X171-C17 X171-C2 X171-C16 | rotatable
- t2: `15 24 14 13`  #X171-C16 X171-N1 X171-C15 X171-C14 | rotatable
- t3: `10 31 11 18`  #X171-C11 X171-S1 X171-C12 X171-C19 | rotatable
- t4: `11 31 10 9`  #X171-C12 X171-S1 X171-C11 X171-C10 | rotatable
- t5: `12 0 28 49`  #X171-C13 X171-C1 X171-O2 X171-H18 | polar-H O-H
- t6: `12 13 29 50`  #X171-C13 X171-C14 X171-O3 X171-H19 | polar-H O-H
- t7: `11 18 30 51`  #X171-C12 X171-C19 X171-O4 X171-H20 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597279, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 349 clusters | 3.6, 3.58, 2.57, 2.42, 2.35, 1.61, 1.54, 1.39, 1.37, 1.36 ... (+339 more, 78.3 % together) |
| 6QLP | 8 clusters | 49.64, 32.28, 17.27, 0.19, 0.19, 0.18, 0.14, 0.11 |
