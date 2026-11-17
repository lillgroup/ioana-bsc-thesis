# 6qmk -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224183, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/6qmk-water` (REST2-run1) | `6qmk.input.sdf` |
| 6QMK complex | `REST2-rotbonds2/6qmk-complex` | `6qmk.input.sdf` |

Torsions t0..t6 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `3 2 1 0`  #X301-C3 X301-C2 X301-O1 X301-C1 | rotatable
- t1: `13 12 4 3`  #X301-C10 X301-C9 X301-C4 X301-C3 | rotatable
- t2: `14 13 12 4`  #X301-C11 X301-C10 X301-C9 X301-C4 | rotatable
- t3: `15 14 13 12`  #X301-O2 X301-C11 X301-C10 X301-C9 | rotatable
- t4: `18 17 12 4`  #X301-C13 X301-C12 X301-C9 X301-C4 | rotatable
- t5: `24 23 22 20`  #X301-N4 X301-C18 X301-C17 X301-C15 | rotatable
- t6: `25 24 23 22`  #X301-C19 X301-N4 X301-C18 X301-C17 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224184, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 97 clusters | 9.71, 8.42, 6.49, 6.05, 6.02, 5.87, 5.7, 4.53, 3.89, 3.6 ... (+87 more, 39.7 % together) |
| 6QMK | 3 clusters | 99.08, 0.7, 0.21 |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224185, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 52 clusters | 9.71, 8.76, 6.49, 6.16, 6.02, 5.87, 5.7, 4.53, 4.14, 3.89 ... (+42 more, 38.7 % together) |
| 6QMK | 2 clusters | 99.08, 0.92 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224186, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 38 clusters | 9.71, 8.86, 7.34, 6.66, 6.22, 6.16, 5.7, 4.53, 4.14, 4.0 ... (+28 more, 36.7 % together) |
| 6QMK | 2 clusters | 99.08, 0.92 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
