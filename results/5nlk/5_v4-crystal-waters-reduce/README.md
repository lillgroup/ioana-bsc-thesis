# 5nlk -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224223, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/5nlk-water` (REST2-run1) | `5nlk.input.sdf` |
| 5NLK complex | `REST2-rotbonds2/5nlk-complex` | `5nlk.input.sdf` |

Torsions t0..t10 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `6 19 14 2`  #X361-C7 X361-C20 X361-C15 X361-C3 | rotatable
- t1: `9 20 15 26`  #X361-C10 X361-C21 X361-C16 X361-N3 | rotatable
- t2: `5 23 16 25`  #X361-C6 X361-C24 X361-C17 X361-N2 | rotatable
- t3: `13 24 17 8`  #X361-C14 X361-N1 X361-C18 X361-C9 | rotatable
- t4: `16 25 18 8`  #X361-C17 X361-N2 X361-C19 X361-C9 | rotatable
- t5: `15 26 21 11`  #X361-C16 X361-N3 X361-C22 X361-C12 | rotatable
- t6: `12 31 22 7`  #X361-C13 X361-O5 X361-C23 X361-C8 | rotatable
- t7: `22 31 12 0`  #X361-C23 X361-O5 X361-C13 X361-C1 | rotatable
- t8: `17 24 13 1`  #X361-C18 X361-N1 X361-C14 X361-C2 | rigid amide-like C-N
- t9: `18 25 16 23`  #X361-C19 X361-N2 X361-C17 X361-C24 | rigid amide-like C-N
- t10: `21 26 15 20`  #X361-C22 X361-N3 X361-C16 X361-C21 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224224, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 254 clusters | 11.15, 10.3, 5.94, 5.24, 5.06, 4.94, 3.36, 3.02, 1.88, 1.8 ... (+244 more, 47.3 % together) |
| 5NLK | 46 clusters | 28.49, 23.3, 13.53, 7.89, 4.96, 3.24, 3.08, 2.62, 2.27, 1.57 ... (+36 more, 9.0 % together) |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224225, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 105 clusters | 13.03, 10.3, 6.75, 5.24, 5.06, 4.94, 3.68, 3.36, 3.02, 2.57 ... (+95 more, 42.1 % together) |
| 5NLK | 31 clusters | 31.12, 25.91, 13.53, 8.38, 4.96, 3.24, 3.08, 1.57, 1.43, 1.35 ... (+21 more, 5.4 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224226, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 78 clusters | 14.77, 13.72, 7.63, 5.81, 5.64, 5.24, 4.24, 3.36, 3.02, 2.97 ... (+68 more, 33.6 % together) |
| 5NLK | 28 clusters | 36.08, 28.98, 14.57, 8.38, 3.24, 1.57, 1.43, 1.35, 0.84, 0.46 ... (+18 more, 3.1 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
