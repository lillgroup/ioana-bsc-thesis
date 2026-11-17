# 4luz -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224215, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4luz-water` (REST2-run1) | `4luz.input.sdf` |
| 4LUZ complex | `REST2-rotbonds2/4luz-complex` | `4luz.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `4 3 0 1`  #X351-C3 X351-C2 X351-C1 X351-O1 | rotatable
- t1: `9 8 6 5`  #X351-C7 X351-C6 X351-C5 X351-C4 | rotatable
- t2: `15 14 11 10`  #X351-O4 X351-C12 X351-C9 X351-C8 | rotatable
- t3: `16 15 14 11`  #X351-C13 X351-O4 X351-C12 X351-C9 | rotatable
- t4: `17 16 15 14`  #X351-C14 X351-C13 X351-O4 X351-C12 | rotatable
- t5: `23 22 19 18`  #X351-C20 X351-C19 X351-C16 X351-C15 | rotatable
- t6: `28 27 26 22`  #X351-C23 X351-C22 X351-N2 X351-C19 | rotatable
- t7: `35 34 24 23`  #X351-O5 X351-C29 X351-C21 X351-C20 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224216, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 98 clusters | 6.32, 6.14, 5.99, 5.98, 5.88, 5.87, 5.77, 5.76, 5.47, 5.38 ... (+88 more, 41.5 % together) |
| 4LUZ | 31 clusters | 86.27, 4.54, 3.12, 2.57, 0.21, 0.19, 0.17, 0.17, 0.15, 0.15 ... (+21 more, 2.4 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224217, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 71 clusters | 6.56, 6.23, 6.18, 6.14, 5.88, 5.87, 5.77, 5.76, 5.66, 5.47 ... (+61 more, 40.5 % together) |
| 4LUZ | 24 clusters | 86.27, 4.54, 3.12, 2.57, 0.57, 0.33, 0.31, 0.21, 0.21, 0.19 ... (+14 more, 1.7 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224218, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 10 clusters | 27.01, 26.07, 22.97, 21.29, 0.76, 0.56, 0.47, 0.32, 0.29, 0.26 |
| 4LUZ | 14 clusters | 86.38, 4.63, 3.39, 3.24, 0.78, 0.33, 0.31, 0.19, 0.15, 0.14 ... (+4 more, 0.4 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
