# 5mg2 -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224134, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/5mg2-water` (REST2-run1) | `5mg2.input.sdf` |
| 5MG2 complex | `REST2-rotbonds2/5mg2-complex` | `5mg2.input.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `24 5 4 23`  #X211-O3 X211-C6 X211-C5 X211-C21 | rotatable
- t1: `9 20 10 11`  #X211-C10 X211-N1 X211-C11 X211-C12 | rotatable
- t2: `22 23 4 5`  #X211-C20 X211-C21 X211-C5 X211-C6 | rotatable
- t3: `4 23 22 8`  #X211-C5 X211-C21 X211-C20 X211-C9 | rotatable
- t4: `4 5 24 47`  #X211-C5 X211-C6 X211-O3 X211-H16 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224135, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 57 clusters | 7.72, 7.33, 6.71, 6.09, 5.19, 5.01, 4.83, 4.49, 3.71, 3.48 ... (+47 more, 45.4 % together) |
| 5MG2 | 56 clusters | 12.4, 7.31, 7.13, 7.02, 6.25, 4.89, 4.54, 4.46, 3.99, 3.72 ... (+46 more, 38.3 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224136, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 50 clusters | 7.72, 7.33, 6.71, 6.2, 5.19, 5.01, 4.93, 4.64, 3.71, 3.48 ... (+40 more, 45.1 % together) |
| 5MG2 | 53 clusters | 12.4, 7.31, 7.13, 7.02, 6.25, 4.89, 4.54, 4.46, 3.99, 3.72 ... (+43 more, 38.3 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224137, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 49 clusters | 7.72, 7.33, 6.71, 6.2, 5.19, 5.01, 4.93, 4.64, 3.71, 3.48 ... (+39 more, 45.1 % together) |
| 5MG2 | 51 clusters | 12.4, 7.31, 7.13, 7.02, 6.58, 4.89, 4.54, 4.46, 3.99, 3.72 ... (+41 more, 37.9 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
