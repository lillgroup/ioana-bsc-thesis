# 4erf-reduce -- 5_v4-reduce

Made by `analysis/report/report.sbatch` (Slurm job 31533375, 2026-11-12, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4erf-water` (REST2-run1) | `4erf.input.sdf` |
| 4ERF complex | `REST2-rotbonds-reduce/4erf-complex` | `4erf.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `8 7 6 0`  #X111-C8 X111-C7 X111-C6 X111-C1 | rotatable
- t1: `8 9 10 11`  #X111-C8 X111-C9 X111-C10 X111-C11 | rotatable
- t2: `12 11 10 9`  #X111-O1 X111-C11 X111-C10 X111-C9 | rotatable
- t3: `19 18 17 15`  #X111-C15 X111-C14 X111-N1 X111-C13 | rotatable
- t4: `20 19 18 17`  #X111-C16 X111-C15 X111-C14 X111-N1 | rotatable
- t5: `22 21 18 17`  #X111-O4 X111-C17 X111-C14 X111-N1 | rotatable
- t6: `26 25 24 7`  #X111-C21 X111-C20 X111-C19 X111-C7 | rotatable
- t7: `18 21 22 51`  #X111-C14 X111-C17 X111-O4 X111-H20 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 31533376, 2026-11-12);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 49 clusters | 41.09, 27.04, 5.39, 4.93, 3.48, 3.2, 2.81, 1.57, 1.38, 1.08 ... (+39 more, 8.0 % together) |
| 4ERF | 11 clusters | 42.88, 33.98, 10.95, 4.86, 3.62, 1.05, 1.03, 0.87, 0.41, 0.19 ... (+1 more, 0.1 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 31533377, 2026-11-12);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 42 clusters | 41.27, 27.04, 5.39, 5.24, 3.48, 3.2, 2.81, 1.57, 1.38, 1.08 ... (+32 more, 7.5 % together) |
| 4ERF | 9 clusters | 43.91, 33.98, 11.1, 4.86, 3.62, 1.05, 0.87, 0.41, 0.19 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 31533378, 2026-11-12);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 39 clusters | 46.51, 30.24, 5.39, 3.48, 2.81, 1.57, 1.38, 1.08, 0.78, 0.74 ... (+29 more, 6.0 % together) |
| 4ERF | 9 clusters | 43.91, 33.98, 11.1, 4.86, 3.62, 1.05, 0.87, 0.41, 0.19 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
