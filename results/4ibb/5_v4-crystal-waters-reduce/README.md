# 4ibb -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224167, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4ibb-water` (REST2-run1) | `4ibb.input.sdf` |
| 4IBB complex | `REST2-rotbonds2/4ibb-complex` | `4ibb.input.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `10 9 8 1`  #X271-C11 X271-C10 X271-C9 X271-C2 | rotatable
- t1: `12 11 10 9`  #X271-O1 X271-C12 X271-C11 X271-C10 | rotatable
- t2: `14 13 11 10`  #X271-C14 X271-C13 X271-C12 X271-C11 | rotatable
- t3: `24 23 22 9`  #X271-C20 X271-C19 X271-N1 X271-C10 | rotatable
- t4: `30 29 28 25`  #X271-C26 X271-C25 X271-C24 X271-C21 | rotatable
- t5: `31 30 29 28`  #X271-O4 X271-C26 X271-C25 X271-C24 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224168, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 29 clusters | 27.23, 18.21, 16.5, 15.72, 8.53, 6.53, 2.2, 0.88, 0.53, 0.4 ... (+19 more, 3.3 % together) |
| 4IBB | 23 clusters | 53.21, 23.09, 9.68, 6.67, 2.96, 0.98, 0.91, 0.3, 0.24, 0.18 ... (+13 more, 1.8 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224169, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 8 clusters | 27.23, 18.21, 16.5, 15.72, 8.53, 6.53, 3.68, 3.61 |
| 4IBB | 8 clusters | 60.96, 23.09, 9.68, 3.26, 1.18, 0.98, 0.74, 0.11 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224170, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 4 clusters | 45.44, 32.22, 15.05, 7.29 |
| 4IBB | 7 clusters | 60.96, 32.77, 3.26, 1.18, 0.98, 0.74, 0.11 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
