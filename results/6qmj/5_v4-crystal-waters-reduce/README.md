# 6qmj -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224269, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/6qmj-water` (REST2-run1) | `6qmj.input.sdf` |
| 6QMJ complex | `REST2-rotbonds2/6qmj-complex` | `6qmj.input.sdf` |

Torsions t0..t8 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `8 7 2 3`  #X411-C9 X411-C8 X411-C3 X411-C4 | rotatable
- t1: `9 8 7 2`  #X411-C10 X411-C9 X411-C8 X411-C3 | rotatable
- t2: `10 9 8 7`  #X411-O1 X411-C10 X411-C9 X411-C8 | rotatable
- t3: `28 15 14 12`  #X411-N4 X411-C15 X411-C14 X411-C12 | rotatable
- t4: `0 20 1 6`  #X411-C1 X411-O2 X411-C2 X411-C7 | rotatable
- t5: `11 26 7 2`  #X411-C11 X411-C21 X411-C8 X411-C3 | rotatable
- t6: `16 28 15 14`  #X411-C16 X411-N4 X411-C15 X411-C14 | rotatable
- t7: `30 29 28 15`  #X411-O4 X411-S1 X411-N4 X411-C15 | rotatable
- t8: `17 32 29 28`  #X411-C17 X411-C23 X411-S1 X411-N4 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224270, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 294 clusters | 4.67, 4.52, 4.45, 3.31, 3.09, 3.02, 2.82, 2.78, 2.69, 2.59 ... (+284 more, 66.0 % together) |
| 6QMJ | 2 clusters | 99.88, 0.12 |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224271, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 60 clusters | 7.52, 7.28, 5.49, 5.35, 4.86, 4.75, 4.72, 4.52, 4.37, 3.86 ... (+50 more, 47.3 % together) |
| 6QMJ | 2 clusters | 99.88, 0.12 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224272, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 21 clusters | 14.83, 13.49, 7.41, 7.28, 5.97, 5.68, 5.67, 5.35, 4.77, 4.75 ... (+11 more, 24.8 % together) |
| 6QMJ | 2 clusters | 99.88, 0.12 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
