# 4x5p -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224151, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4x5p-water` (REST2-run1) | `4x5p.input.sdf` |
| 4X5P complex | `REST2-rotbonds2/4x5p-complex` | `4x5p.input.sdf` |

Torsions t0..t10 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `6 5 4 19`  #X251-C5 X251-C4 X251-C3 X251-O8 | rotatable
- t1: `0 13 9 7`  #X251-O1 X251-C10 X251-C7 X251-C6 | rotatable
- t2: `1 21 20 4`  #X251-C1 X251-C13 X251-N1 X251-C3 | rotatable
- t3: `15 26 25 2`  #X251-C11 X251-O9 X251-C16 X251-C2 | rotatable
- t4: `25 26 15 3`  #X251-C16 X251-O9 X251-C11 X251-O2 | rotatable
- t5: `16 30 17 3`  #X251-O6 X251-C20 X251-C12 X251-O2 | rotatable
- t6: `21 20 4 5`  #X251-C13 X251-N1 X251-C3 X251-C4 | rigid amide-like C-N
- t7: `15 27 14 38`  #X251-C11 X251-C17 X251-O5 X251-H8 | polar-H O-H
- t8: `27 28 10 35`  #X251-C17 X251-C18 X251-O4 X251-H5 | polar-H O-H
- t9: `17 29 8 34`  #X251-C12 X251-C19 X251-O3 X251-H4 | polar-H O-H
- t10: `17 30 16 40`  #X251-C12 X251-C20 X251-O6 X251-H10 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224152, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 661 clusters | 0.63, 0.61, 0.56, 0.56, 0.5, 0.47, 0.4, 0.38, 0.37, 0.35 ... (+651 more, 95.2 % together) |
| 4X5P | 18 clusters | 22.85, 18.93, 18.62, 18.53, 4.9, 4.45, 4.26, 4.09, 0.78, 0.69 ... (+8 more, 1.9 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224153, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 533 clusters | 0.79, 0.63, 0.63, 0.61, 0.6, 0.58, 0.56, 0.56, 0.56, 0.56 ... (+523 more, 93.9 % together) |
| 4X5P | 10 clusters | 27.3, 23.43, 23.2, 22.71, 1.0, 0.78, 0.69, 0.66, 0.12, 0.12 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224154, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 373 clusters | 1.08, 0.89, 0.87, 0.81, 0.77, 0.75, 0.65, 0.64, 0.63, 0.61 ... (+363 more, 92.3 % together) |
| 4X5P | 9 clusters | 27.3, 23.43, 23.2, 22.71, 1.0, 0.81, 0.78, 0.66, 0.12 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
