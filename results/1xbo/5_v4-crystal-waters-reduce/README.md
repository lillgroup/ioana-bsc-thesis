# 1xbo -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224231, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/1xbo-water` (REST2-run1) | `1xbo.input.sdf` |
| 1XBO complex | `REST2-rotbonds2/1xbo-complex` | `1xbo.input.sdf` |

Torsions t0..t11 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `7 6 1 0`  #X381-O1 X381-C7 X381-C2 X381-C1 | rotatable
- t1: `12 11 9 8`  #X381-O2 X381-C10 X381-C8 X381-N1 | rotatable
- t2: `28 14 3 2`  #X381-C21 X381-C11 X381-C4 X381-C3 | rotatable
- t3: `22 21 20 15`  #X381-O4 X381-C18 X381-C17 X381-C12 | rotatable
- t4: `27 26 15 16`  #X381-C20 X381-O7 X381-C12 X381-C13 | rotatable
- t5: `28 27 26 15`  #X381-C21 X381-C20 X381-O7 X381-C12 | rotatable
- t6: `14 28 27 26`  #X381-C11 X381-C21 X381-C20 X381-O7 | rotatable
- t7: `30 29 10 6`  #X381-O8 X381-C22 X381-C9 X381-C7 | rotatable
- t8: `24 23 21 20`  #X381-C19 X381-O5 X381-C18 X381-C17 | rigid ester C(=O)-O
- t9: `27 28 14 3`  #X381-C20 X381-C21 X381-C11 X381-C4 | rigid double bond
- t10: `18 19 25 42`  #X381-C15 X381-C16 X381-O6 X381-H12 | polar-H O-H
- t11: `10 29 30 48`  #X381-C9 X381-C22 X381-O8 X381-H18 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224246, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 671 clusters | 0.46, 0.39, 0.39, 0.38, 0.36, 0.35, 0.33, 0.32, 0.32, 0.31 ... (+661 more, 96.2 % together) |
| 1XBO | 25 clusters | 46.13, 30.84, 10.37, 5.17, 2.79, 0.93, 0.5, 0.31, 0.3, 0.26 ... (+15 more, 2.4 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224247, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 491 clusters | 0.99, 0.95, 0.89, 0.83, 0.81, 0.8, 0.77, 0.72, 0.69, 0.66 ... (+481 more, 91.8 % together) |
| 1XBO | 17 clusters | 76.98, 10.87, 7.96, 1.08, 0.65, 0.42, 0.31, 0.3, 0.26, 0.18 ... (+7 more, 1.0 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224248, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 373 clusters | 3.53, 1.42, 1.05, 1.0, 0.99, 0.98, 0.9, 0.89, 0.89, 0.87 ... (+363 more, 87.4 % together) |
| 1XBO | 12 clusters | 88.0, 7.96, 2.06, 0.42, 0.31, 0.3, 0.26, 0.18, 0.17, 0.13 ... (+2 more, 0.2 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
