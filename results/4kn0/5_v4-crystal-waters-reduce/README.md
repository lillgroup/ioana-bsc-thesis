# 4kn0 -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224277, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4kn0-water` (REST2-run1) | `4kn0.input.sdf` |
| 4KN0 complex | `REST2-rotbonds2/4kn0-complex` | `4kn0.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `13 12 8 7`  #X421-N7 X421-C7 X421-C4 X421-N5 | rotatable
- t1: `14 13 12 8`  #X421-C8 X421-N7 X421-C7 X421-C4 | rotatable
- t2: `17 18 13 12`  #X421-C11 X421-C12 X421-N7 X421-C7 | rotatable
- t3: `22 21 15 16`  #X421-O1 X421-C15 X421-C9 X421-C10 | rotatable
- t4: `25 24 23 21`  #X421-C17 X421-C16 X421-N8 X421-C15 | rotatable
- t5: `26 25 24 23`  #X421-O2 X421-C17 X421-C16 X421-N8 | rotatable
- t6: `29 28 24 23`  #X421-C19 X421-C18 X421-C16 X421-N8 | rotatable
- t7: `30 29 28 24`  #X421-C20 X421-C19 X421-C18 X421-C16 | rotatable
- t8: `31 30 29 28`  #X421-O4 X421-C20 X421-C19 X421-C18 | rotatable
- t9: `24 23 21 15`  #X421-C16 X421-N8 X421-C15 X421-C9 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224278, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 320 clusters | 10.91, 10.7, 7.35, 6.28, 6.03, 4.65, 2.47, 2.29, 2.04, 1.43 ... (+310 more, 45.8 % together) |
| 4KN0 | 41 clusters | 50.26, 20.14, 11.04, 4.26, 3.98, 2.35, 1.84, 1.18, 0.6, 0.21 ... (+31 more, 4.1 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224279, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 97 clusters | 17.19, 16.72, 11.17, 7.04, 5.36, 3.7, 2.63, 2.44, 2.37, 2.29 ... (+87 more, 29.1 % together) |
| 4KN0 | 20 clusters | 50.88, 20.14, 11.04, 6.83, 4.35, 1.96, 1.79, 0.92, 0.25, 0.25 ... (+10 more, 1.6 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224280, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 67 clusters | 17.19, 16.72, 11.17, 7.04, 5.36, 3.7, 2.94, 2.63, 2.44, 2.29 ... (+57 more, 28.5 % together) |
| 4KN0 | 15 clusters | 58.23, 20.14, 11.04, 4.46, 1.96, 1.79, 0.92, 0.25, 0.25, 0.21 ... (+5 more, 0.8 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
