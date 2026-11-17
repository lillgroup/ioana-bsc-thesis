# 5fnt -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33224199, 2026-11-16, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/5fnt-water` (REST2-run1) | `5fnt.input.sdf` |
| 5FNT complex | `REST2-rotbonds2/5fnt-complex` | `5fnt.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `3 2 1 0`  #X321-C3 X321-C2 X321-N1 X321-C1 | rotatable
- t1: `4 3 2 1`  #X321-C4 X321-C3 X321-C2 X321-N1 | rotatable
- t2: `11 10 5 4`  #X321-C10 X321-C9 X321-C5 X321-C4 | rotatable
- t3: `12 11 10 5`  #X321-C11 X321-C10 X321-C9 X321-C5 | rotatable
- t4: `13 12 11 10`  #X321-O1 X321-C11 X321-C10 X321-C9 | rotatable
- t5: `16 15 10 5`  #X321-C13 X321-C12 X321-C9 X321-C5 | rotatable
- t6: `26 25 1 0`  #X321-O3 X321-S1 X321-N1 X321-C1 | rotatable
- t7: `29 28 25 1`  #X321-C20 X321-C19 X321-S1 X321-N1 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224200, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 181 clusters | 6.01, 5.8, 4.78, 4.24, 3.86, 3.52, 3.4, 3.34, 3.08, 3.06 ... (+171 more, 58.9 % together) |
| 5FNT | 3 clusters | 99.23, 0.64, 0.13 |

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224201, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 60 clusters | 6.18, 6.07, 5.39, 4.61, 4.52, 4.37, 3.86, 3.76, 3.69, 3.51 ... (+50 more, 54.0 % together) |
| 5FNT | 2 clusters | 99.23, 0.77 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33224202, 2026-11-16);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 37 clusters | 6.18, 6.07, 5.39, 4.71, 4.52, 4.37, 3.99, 3.9, 3.8, 3.76 ... (+27 more, 53.3 % together) |
| 5FNT | 2 clusters | 99.23, 0.77 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
