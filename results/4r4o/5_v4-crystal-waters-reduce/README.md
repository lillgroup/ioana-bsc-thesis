# 4r4o -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224207, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4r4o-water` (REST2-run1) | `4r4o.input.sdf` |
| 4R4O complex | `REST2-rotbonds2/4r4o-complex` | `4r4o.input.sdf` |

Torsions t0..t8 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `25 5 4 24`  #X341-C20 X341-C6 X341-C5 X341-N3 | rotatable
- t1: `12 11 10 9`  #X341-N1 X341-C11 X341-C10 X341-O1 | rotatable
- t2: `1 14 13 12`  #X341-C2 X341-C12 X341-N2 X341-N1 | rotatable
- t3: `13 20 2 21`  #X341-N2 X341-C16 X341-C3 X341-C17 | rotatable
- t4: `24 23 22 3`  #X341-N3 X341-C19 X341-C18 X341-C4 | rotatable
- t5: `23 24 4 5`  #X341-C19 X341-N3 X341-C5 X341-C6 | rotatable
- t6: `6 26 7 27`  #X341-C7 X341-C21 X341-C8 X341-C22 | rotatable
- t7: `30 29 28 8`  #X341-O2 X341-C24 X341-C23 X341-C9 | rotatable
- t8: `4 24 23 22`  #X341-C5 X341-N3 X341-C19 X341-C18 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224208, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 159 clusters | 3.6, 3.31, 3.16, 3.01, 3.01, 2.99, 2.96, 2.89, 2.76, 2.71 ... (+149 more, 69.7 % together) |
| 4R4O | 34 clusters | 66.95, 21.63, 2.25, 1.86, 1.74, 1.2, 0.46, 0.29, 0.21, 0.2 ... (+24 more, 3.2 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224209, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 24 clusters | 8.69, 8.33, 7.38, 7.38, 6.58, 6.2, 6.16, 5.66, 5.23, 5.14 ... (+14 more, 33.3 % together) |
| 4R4O | 23 clusters | 66.95, 21.78, 2.25, 1.86, 1.74, 1.33, 0.86, 0.74, 0.46, 0.24 ... (+13 more, 1.7 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224210, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 14 clusters | 15.71, 10.8, 8.69, 8.38, 7.38, 6.64, 6.58, 6.2, 5.66, 5.41 ... (+4 more, 18.6 % together) |
| 4R4O | 13 clusters | 74.74, 21.78, 0.86, 0.74, 0.46, 0.46, 0.19, 0.17, 0.13, 0.12 ... (+3 more, 0.3 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
