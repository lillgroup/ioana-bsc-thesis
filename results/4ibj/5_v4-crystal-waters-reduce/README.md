# 4ibj -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224143, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/4ibj-water` (REST2-run1) | `4ibj.input.sdf` |
| 4IBJ complex | `REST2-rotbonds2/4ibj-complex` | `4ibj.input.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `11 10 9 5`  #X221-N1 X221-C8 X221-C7 X221-C3 | rotatable
- t1: `13 12 11 10`  #X221-C10 X221-C9 X221-N1 X221-C8 | rotatable
- t2: `16 15 14 13`  #X221-O1 X221-C12 X221-C11 X221-C10 | rotatable
- t3: `27 26 25 10`  #X221-O5 X221-C19 X221-C18 X221-C8 | rotatable
- t4: `29 28 26 25`  #X221-C21 X221-C20 X221-C19 X221-C18 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224144, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 18 clusters | 23.76, 11.21, 10.83, 10.58, 9.18, 6.39, 5.79, 5.0, 4.13, 3.74 ... (+8 more, 9.4 % together) |
| 4IBJ | 1 clusters | 100.0 |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224145, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 16 clusters | 23.76, 11.21, 10.83, 10.58, 9.18, 6.39, 5.79, 5.0, 4.13, 3.74 ... (+6 more, 9.4 % together) |
| 4IBJ | 1 clusters | 100.0 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224146, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 16 clusters | 23.76, 11.21, 10.83, 10.58, 9.18, 6.39, 5.79, 5.0, 4.13, 3.74 ... (+6 more, 9.4 % together) |
| 4IBJ | 1 clusters | 100.0 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
