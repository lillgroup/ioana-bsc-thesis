# 3qxt -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224285, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/3qxt-water` (REST2-run1) | `3qxt.input.sdf` |
| 3QXT complex | `REST2-rotbonds2/3qxt-complex` | `3qxt.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `13 12 8 7`  #X431-N7 X431-C7 X431-C4 X431-N5 | rotatable
- t1: `14 13 12 8`  #X431-C8 X431-N7 X431-C7 X431-C4 | rotatable
- t2: `17 18 13 12`  #X431-C11 X431-C12 X431-N7 X431-C7 | rotatable
- t3: `22 21 15 16`  #X431-O1 X431-C15 X431-C9 X431-C10 | rotatable
- t4: `25 24 23 21`  #X431-C17 X431-C16 X431-N8 X431-C15 | rotatable
- t5: `26 25 24 23`  #X431-O2 X431-C17 X431-C16 X431-N8 | rotatable
- t6: `29 28 24 23`  #X431-C19 X431-C18 X431-C16 X431-N8 | rotatable
- t7: `30 29 28 24`  #X431-C20 X431-C19 X431-C18 X431-C16 | rotatable
- t8: `31 30 29 28`  #X431-O4 X431-C20 X431-C19 X431-C18 | rotatable
- t9: `24 23 21 15`  #X431-C16 X431-N8 X431-C15 X431-C9 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224286, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 271 clusters | 9.39, 9.17, 7.89, 6.43, 5.97, 4.71, 4.32, 4.06, 3.78, 2.95 ... (+261 more, 41.4 % together) |
| 3QXT | 75 clusters | 22.23, 20.29, 17.29, 13.86, 6.97, 4.09, 3.73, 1.42, 0.63, 0.43 ... (+65 more, 9.0 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224287, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 121 clusters | 15.36, 13.23, 9.88, 7.89, 4.71, 4.46, 4.32, 2.84, 2.39, 2.31 ... (+111 more, 32.6 % together) |
| 3QXT | 42 clusters | 22.23, 20.29, 17.29, 13.86, 7.27, 4.46, 3.96, 1.8, 1.5, 1.07 ... (+32 more, 6.2 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224288, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 59 clusters | 15.36, 13.23, 12.21, 9.88, 4.71, 4.46, 2.84, 2.84, 2.79, 2.53 ... (+49 more, 29.1 % together) |
| 3QXT | 34 clusters | 22.23, 20.29, 17.29, 13.86, 8.14, 4.46, 3.96, 1.93, 1.5, 1.07 ... (+24 more, 5.2 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
