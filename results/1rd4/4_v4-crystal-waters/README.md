# 1rd4 -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620734, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/1rd4-water` (REST2-run1) | `1rd4.input.sdf` |
| 1RD4 complex | `REST2-rotbonds/1rd4-complex` | `1rd4.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `3 2 1 0`  #X131-C2 X131-C1 X131-N1 X131-O1 | rotatable
- t1: `9 8 6 5`  #X131-C8 X131-C7 X131-C5 X131-C4 | rotatable
- t2: `14 13 12 11`  #X131-C11 X131-N3 X131-C10 X131-N2 | rotatable
- t3: `24 23 3 2`  #X131-C18 X131-S1 X131-C2 X131-C1 | rotatable
- t4: `25 24 23 3`  #X131-C19 X131-C18 X131-S1 X131-C2 | rotatable
- t5: `31 30 25 24`  #X131-C24 X131-O3 X131-C19 X131-C18 | rotatable
- t6: `32 31 30 25`  #X131-C25 X131-C24 X131-O3 X131-C19 | rotatable
- t7: `20 19 16 15`  #X131-O2 X131-C15 X131-N4 X131-C12 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597290, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 304 clusters | 2.84, 2.74, 2.7, 2.63, 2.62, 2.41, 2.19, 2.12, 2.11, 2.08 ... (+294 more, 75.5 % together) |
| 1RD4 | 68 clusters | 33.95, 19.57, 8.6, 7.02, 5.62, 4.83, 3.46, 1.84, 1.22, 1.1 ... (+58 more, 12.8 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857840, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 157 clusters | 6.18, 5.2, 5.14, 5.01, 4.97, 4.96, 4.9, 4.65, 4.63, 3.89 ... (+147 more, 50.4 % together) |
| 1RD4 | 28 clusters | 39.57, 21.41, 8.8, 7.02, 5.17, 3.94, 2.31, 1.53, 1.26, 1.23 ... (+18 more, 7.8 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857841, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 40 clusters | 8.25, 8.08, 7.88, 7.7, 7.23, 7.14, 5.95, 5.8, 5.59, 5.0 ... (+30 more, 31.4 % together) |
| 1RD4 | 8 clusters | 67.52, 23.18, 3.94, 1.53, 1.22, 1.16, 0.9, 0.56 |
