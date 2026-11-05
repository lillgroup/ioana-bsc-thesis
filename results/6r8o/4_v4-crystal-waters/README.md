# 6r8o -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620729, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/6r8o-water` (REST2-run1) | `6r8o.input.sdf` |
| 6R8O complex | `REST2-rotbonds/6r8o-complex` | `6r8o.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 20 1 2`  #X081-C1 X081-S1 X081-C2 X081-C3 | rotatable
- t1: `9 21 6 1`  #X081-N1 X081-C20 X081-C7 X081-C2 | rotatable
- t2: `9 23 10 25`  #X081-N1 X081-C22 X081-C10 X081-N2 | rotatable
- t3: `11 25 10 23`  #X081-C11 X081-N2 X081-C10 X081-C22 | rotatable
- t4: `11 27 12 28`  #X081-C11 X081-N3 X081-C12 X081-C23 | rotatable
- t5: `13 28 12 27`  #X081-C13 X081-C23 X081-C12 X081-N3 | rotatable
- t6: `10 23 9 8`  #X081-C10 X081-C22 X081-N1 X081-C9 | rigid amide-like C-N
- t7: `10 25 11 26`  #X081-C10 X081-N2 X081-C11 X081-O2 | rigid amide-like C-N
- t8: `12 27 11 25`  #X081-C12 X081-N3 X081-C11 X081-N2 | rigid amide-like C-N
- t9: `15 18 33 66`  #X081-C15 X081-C18 X081-O3 X081-H32 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597280, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 413 clusters | 1.82, 1.52, 1.45, 1.44, 1.41, 1.38, 1.2, 1.17, 1.14, 1.03 ... (+403 more, 86.4 % together) |
| 6R8O | 12 clusters | 51.58, 15.15, 14.24, 7.13, 4.21, 2.97, 2.76, 0.81, 0.59, 0.32 ... (+2 more, 0.2 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
