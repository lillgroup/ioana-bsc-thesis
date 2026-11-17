# 4jfm -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224261, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4jfm-water` (REST2-run1) | `4jfm.input.sdf` |
| 4JFM complex | `REST2-rotbonds2/4jfm-complex` | `4jfm.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `3 2 1 0`  #X401-C3 X401-C2 X401-O1 X401-C1 | rotatable
- t1: `6 5 4 2`  #X401-C5 X401-O2 X401-C4 X401-C2 | rotatable
- t2: `11 10 9 3`  #X401-C9 X401-O3 X401-C8 X401-C3 | rotatable
- t3: `12 11 10 9`  #X401-C10 X401-C9 X401-O3 X401-C8 | rotatable
- t4: `13 12 11 10`  #X401-O4 X401-C10 X401-C9 X401-O3 | rotatable
- t5: `14 13 12 11`  #X401-C11 X401-O4 X401-C10 X401-C9 | rotatable
- t6: `17 16 14 13`  #X401-C13 X401-C12 X401-C11 X401-O4 | rotatable
- t7: `23 22 21 16`  #X401-O6 X401-S1 X401-N1 X401-C12 | rotatable
- t8: `26 25 22 21`  #X401-C18 X401-C17 X401-S1 X401-N1 | rotatable
- t9: `15 14 13 12`  #X401-O5 X401-C11 X401-O4 X401-C10 | rigid ester C(=O)-O

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224262, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 302 clusters | 5.07, 4.0, 3.55, 3.22, 2.8, 2.46, 2.38, 2.3, 2.28, 2.19 ... (+292 more, 69.7 % together) |
| 4JFM | 134 clusters | 25.72, 11.34, 9.99, 8.03, 4.33, 4.01, 3.95, 3.33, 3.0, 2.11 ... (+124 more, 24.1 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224263, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 155 clusters | 6.41, 4.48, 4.29, 3.63, 3.17, 2.83, 2.8, 2.67, 2.59, 2.56 ... (+145 more, 64.6 % together) |
| 4JFM | 87 clusters | 25.83, 11.45, 10.17, 8.25, 4.8, 4.58, 4.07, 3.83, 3.11, 2.26 ... (+77 more, 21.6 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224264, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 72 clusters | 11.27, 10.16, 8.69, 6.54, 6.22, 4.99, 4.72, 3.32, 2.83, 2.77 ... (+62 more, 38.5 % together) |
| 4JFM | 44 clusters | 25.83, 11.45, 10.17, 8.36, 5.74, 4.91, 4.33, 3.94, 3.62, 2.41 ... (+34 more, 19.2 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
