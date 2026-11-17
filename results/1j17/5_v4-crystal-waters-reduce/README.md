# 1j17 -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224118, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/1j17-water` (REST2-run1) | `1j17.input.sdf` |
| 1J17 complex | `REST2-rotbonds2/1j17-complex` | `1j17.input.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `7 6 5 2`  #X191-C6 X191-N2 X191-C5 X191-C2 | rotatable
- t1: `13 12 11 8`  #X191-O1 X191-C11 X191-C10 X191-C7 | rotatable
- t2: `21 20 19 16`  #X191-O2 X191-S1 X191-N4 X191-C13 | rotatable
- t3: `24 23 20 19`  #X191-C17 X191-C16 X191-S1 X191-N4 | rotatable
- t4: `15 14 12 11`  #X191-C12 X191-N3 X191-C11 X191-C10 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224119, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 24 clusters | 13.67, 12.28, 11.91, 10.94, 10.74, 7.5, 6.43, 6.32, 4.34, 3.27 ... (+14 more, 12.6 % together) |
| 1J17 | 11 clusters | 29.25, 26.16, 14.39, 11.56, 6.54, 4.85, 3.42, 2.28, 0.71, 0.69 ... (+1 more, 0.1 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224120, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 4 clusters | 27.9, 25.09, 23.99, 23.02 |
| 1J17 | 5 clusters | 31.86, 29.25, 16.41, 15.11, 7.36 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224121, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 4 clusters | 27.9, 25.09, 23.99, 23.02 |
| 1J17 | 2 clusters | 54.33, 45.67 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
