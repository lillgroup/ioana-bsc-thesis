# 4jfk -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224191, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4jfk-water` (REST2-run1) | `4jfk.input.sdf` |
| 4JFK complex | `REST2-rotbonds2/4jfk-complex` | `4jfk.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `35 20 17 25`  #X311-N3 X311-C16 X311-C13 X311-O7 | rotatable
- t1: `27 25 17 20`  #X311-C18 X311-O7 X311-C13 X311-C16 | rotatable
- t2: `10 27 25 17`  #X311-C6 X311-C18 X311-O7 X311-C13 | rotatable
- t3: `11 29 23 5`  #X311-C7 X311-C20 X311-O5 X311-C4 | rotatable
- t4: `14 30 24 6`  #X311-C10 X311-C21 X311-O6 X311-C5 | rotatable
- t5: `0 35 20 17`  #X311-C1 X311-N3 X311-C16 X311-C13 | rotatable
- t6: `1 36 28 12`  #X311-N1 X311-S2 X311-C19 X311-C8 | rotatable
- t7: `8 36 1 3`  #X311-O3 X311-S2 X311-N1 X311-C2 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224192, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 158 clusters | 7.83, 6.99, 6.91, 6.8, 6.16, 5.44, 5.02, 4.0, 3.41, 3.37 ... (+148 more, 44.0 % together) |
| 4JFK | 32 clusters | 68.49, 16.46, 8.01, 1.12, 0.89, 0.66, 0.4, 0.34, 0.32, 0.24 ... (+22 more, 3.0 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224193, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 63 clusters | 12.03, 11.33, 8.29, 7.1, 6.94, 6.19, 5.53, 4.24, 3.89, 3.4 ... (+53 more, 31.0 % together) |
| 4JFK | 23 clusters | 68.49, 17.07, 8.01, 1.12, 0.89, 0.78, 0.61, 0.53, 0.45, 0.32 ... (+13 more, 1.7 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224194, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 38 clusters | 23.61, 15.41, 11.72, 7.44, 7.1, 4.37, 3.51, 3.46, 2.51, 2.43 ... (+28 more, 18.4 % together) |
| 4JFK | 22 clusters | 68.62, 17.07, 8.01, 1.12, 0.89, 0.78, 0.61, 0.53, 0.45, 0.32 ... (+12 more, 1.6 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
