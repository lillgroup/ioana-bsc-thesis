# 5h9r -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224159, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/5h9r-water` (REST2-run1) | `5h9r.input.sdf` |
| 5H9R complex | `REST2-rotbonds2/5h9r-complex` | `5h9r.input.sdf` |

Torsions t0..t12 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `12 4 3 8`  #X261-O3 X261-C4 X261-C3 X261-O1 | rotatable
- t1: `1 6 0 7`  #X261-C1 X261-C6 X261-S1 X261-C7 | rotatable
- t2: `9 7 0 6`  #X261-O2 X261-C7 X261-S1 X261-C6 | rotatable
- t3: `13 11 10 9`  #X261-O4 X261-C9 X261-C8 X261-O2 | rotatable
- t4: `24 23 2 1`  #X261-N2 X261-N1 X261-C2 X261-C1 | rotatable
- t5: `29 28 26 25`  #X261-C17 X261-C16 X261-C14 X261-N3 | rotatable
- t6: `3 4 12 44`  #X261-C3 X261-C4 X261-O3 X261-H11 | polar-H O-H
- t7: `10 11 13 45`  #X261-C8 X261-C9 X261-O4 X261-H12 | polar-H O-H
- t8: `2 14 16 48`  #X261-C2 X261-C10 X261-O5 X261-H15 | polar-H O-H
- t9: `10 15 17 49`  #X261-C8 X261-C11 X261-O6 X261-H16 | polar-H O-H
- t10: `15 18 19 51`  #X261-C11 X261-C12 X261-O7 X261-H18 | polar-H O-H
- t11: `2 1 21 53`  #X261-C2 X261-C1 X261-O8 X261-H20 | polar-H O-H
- t12: `7 20 22 54`  #X261-C7 X261-C13 X261-O9 X261-H21 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224160, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 734 clusters | 0.27, 0.23, 0.23, 0.22, 0.22, 0.21, 0.2, 0.2, 0.2, 0.2 ... (+724 more, 97.8 % together) |
| 5H9R | 234 clusters | 3.52, 2.74, 2.69, 2.61, 2.58, 2.54, 2.31, 2.3, 1.99, 1.94 ... (+224 more, 74.7 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224161, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 726 clusters | 0.27, 0.24, 0.23, 0.23, 0.23, 0.22, 0.22, 0.22, 0.21, 0.21 ... (+716 more, 97.6 % together) |
| 5H9R | 201 clusters | 3.52, 2.74, 2.69, 2.61, 2.58, 2.54, 2.31, 2.3, 1.99, 1.94 ... (+191 more, 74.8 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224162, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 590 clusters | 0.36, 0.35, 0.35, 0.34, 0.34, 0.34, 0.34, 0.33, 0.33, 0.33 ... (+580 more, 96.6 % together) |
| 5H9R | 175 clusters | 3.52, 3.14, 2.74, 2.69, 2.58, 2.54, 2.31, 2.3, 1.99, 1.97 ... (+165 more, 74.2 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
