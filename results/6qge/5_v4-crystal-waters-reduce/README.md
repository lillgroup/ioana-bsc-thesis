# 6qge -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224293, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/6qge-water` (REST2-run1) | `6qge.input.sdf` |
| 6QGE complex | `REST2-rotbonds2/6qge-complex` | `6qge.input.sdf` |

Torsions t0..t12 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `15 14 9 8`  #X451-C16 X451-C15 X451-C10 X451-C9 | rotatable
- t1: `31 19 16 0`  #X451-O3 X451-C20 X451-C17 X451-C1 | rotatable
- t2: `2 20 3 26`  #X451-C3 X451-C21 X451-C4 X451-N4 | rotatable
- t3: `3 20 2 34`  #X451-C4 X451-C21 X451-C3 X451-S1 | rotatable
- t4: `15 25 17 0`  #X451-C16 X451-N3 X451-C18 X451-C1 | rotatable
- t5: `21 26 3 20`  #X451-C22 X451-N4 X451-C4 X451-C21 | rotatable
- t6: `1 34 2 20`  #X451-C2 X451-S1 X451-C3 X451-C21 | rotatable
- t7: `2 34 1 18`  #X451-C3 X451-S1 X451-C2 X451-C19 | rotatable
- t8: `21 35 4 5`  #X451-C22 X451-C23 X451-C5 X451-C6 | rotatable
- t9: `16 0 30 59`  #X451-C17 X451-C1 X451-O2 X451-H21 | polar-H O-H
- t10: `16 19 31 60`  #X451-C17 X451-C20 X451-O3 X451-H22 | polar-H O-H
- t11: `1 18 32 61`  #X451-C2 X451-C19 X451-O4 X451-H23 | polar-H O-H
- t12: `2 20 33 62`  #X451-C3 X451-C21 X451-O5 X451-H24 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224294, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 738 clusters | 0.21, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 ... (+728 more, 98.0 % together) |
| 6QGE | 86 clusters | 11.92, 11.43, 9.38, 9.17, 6.39, 4.75, 4.67, 4.17, 4.03, 3.45 ... (+76 more, 30.6 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224295, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 733 clusters | 0.25, 0.22, 0.22, 0.21, 0.21, 0.21, 0.2, 0.2, 0.2, 0.2 ... (+723 more, 97.9 % together) |
| 6QGE | 65 clusters | 23.34, 18.55, 10.15, 9.66, 7.63, 4.03, 3.45, 2.48, 1.87, 1.81 ... (+55 more, 17.0 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224296, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 733 clusters | 0.25, 0.22, 0.22, 0.21, 0.21, 0.21, 0.2, 0.2, 0.2, 0.2 ... (+723 more, 97.9 % together) |
| 6QGE | 61 clusters | 23.34, 18.55, 10.15, 9.66, 7.63, 4.29, 4.03, 3.45, 1.87, 1.57 ... (+51 more, 15.4 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
