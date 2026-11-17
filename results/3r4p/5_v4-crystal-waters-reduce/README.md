# 3r4p -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224175, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/3r4p-water` (REST2-run1) | `3r4p.input.sdf` |
| 3R4P complex | `REST2-rotbonds2/3r4p-complex` | `3r4p.input.sdf` |

Torsions t0..t8 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `26 8 5 0`  #X281-C14 X281-O2 X281-C2 X281-C1 | rotatable
- t1: `11 24 0 5`  #X281-N3 X281-C12 X281-C1 X281-C2 | rotatable
- t2: `27 26 8 5`  #X281-C15 X281-C14 X281-O2 X281-C2 | rotatable
- t3: `17 27 26 8`  #X281-N6 X281-C15 X281-C14 X281-O2 | rotatable
- t4: `26 27 17 19`  #X281-C14 X281-C15 X281-N6 X281-N7 | rotatable
- t5: `32 31 15 25`  #X281-C20 X281-C19 X281-N5 X281-C13 | rotatable
- t6: `1 32 31 15`  #X281-F1 X281-C20 X281-C19 X281-N5 | rotatable
- t7: `3 25 2 20`  #X281-O1 X281-C13 X281-N1 X281-C8 | rigid amide-like C-N
- t8: `2 25 15 31`  #X281-N1 X281-C13 X281-N5 X281-C19 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224176, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 95 clusters | 16.77, 10.58, 6.58, 6.2, 5.07, 4.46, 3.31, 3.22, 2.57, 2.37 ... (+85 more, 38.8 % together) |
| 3R4P | 12 clusters | 93.43, 1.71, 1.64, 0.8, 0.62, 0.58, 0.33, 0.24, 0.21, 0.18 ... (+2 more, 0.2 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224177, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 70 clusters | 16.77, 15.62, 7.68, 6.58, 6.2, 5.07, 3.31, 2.37, 2.29, 2.24 ... (+60 more, 31.8 % together) |
| 3R4P | 11 clusters | 93.43, 1.96, 1.64, 0.8, 0.62, 0.58, 0.33, 0.21, 0.18, 0.12 ... (+1 more, 0.1 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224178, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 52 clusters | 40.29, 12.77, 9.04, 6.96, 3.31, 2.37, 2.12, 2.09, 1.71, 1.66 ... (+42 more, 17.6 % together) |
| 3R4P | 9 clusters | 93.43, 2.4, 1.96, 0.8, 0.62, 0.33, 0.21, 0.12, 0.12 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
