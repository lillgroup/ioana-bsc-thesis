# 1fkh -- 5_v4-crystal-waters-reduce

Made by `analysis/report/report.sbatch` (Slurm job 33225657, 2026-11-17, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds2/1fkh-water` (REST2-run1) | `1fkh.input.sdf` |
| 1FKH complex | `REST2-rotbonds2/1fkh-complex` | `1fkh.input.sdf` |

Torsions t0..t10 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `2 1 0 7`  #X391-C3 X391-C2 X391-C1 X391-O1 | rotatable
- t1: `12 11 9 6`  #X391-O4 X391-C8 X391-C7 X391-N1 | rotatable
- t2: `14 13 11 9`  #X391-C10 X391-C9 X391-C8 X391-C7 | rotatable
- t3: `15 14 13 11`  #X391-C11 X391-C10 X391-C9 X391-C8 | rotatable
- t4: `19 18 8 0`  #X391-C15 X391-C14 X391-O2 X391-C1 | rotatable
- t5: `20 19 18 8`  #X391-C16 X391-C15 X391-C14 X391-O2 | rotatable
- t6: `21 20 19 18`  #X391-C17 X391-C16 X391-C15 X391-C14 | rotatable
- t7: `22 21 20 19`  #X391-C18 X391-C17 X391-C16 X391-C15 | rotatable
- t8: `28 27 18 8`  #X391-C24 X391-C23 X391-C14 X391-O2 | rotatable
- t9: `18 8 0 1`  #X391-C14 X391-O2 X391-C1 X391-C2 | rigid ester C(=O)-O
- t10: `10 9 6 1`  #X391-O3 X391-C7 X391-N1 X391-C2 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33225658, 2026-11-17);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 1 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 342 clusters | 2.55, 2.34, 2.26, 2.25, 2.1, 1.78, 1.69, 1.58, 1.55, 1.54 ... (+332 more, 80.4 % together) |
| 1FKH | 45 clusters | 28.28, 25.97, 6.98, 5.87, 5.68, 5.45, 3.02, 2.55, 1.79, 1.29 ... (+35 more, 13.1 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33225659, 2026-11-17);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 254 clusters | 3.61, 2.86, 2.55, 2.34, 2.25, 1.97, 1.89, 1.78, 1.69, 1.58 ... (+244 more, 77.5 % together) |
| 1FKH | 31 clusters | 28.28, 27.12, 6.98, 6.49, 6.22, 5.87, 3.27, 3.02, 2.06, 1.83 ... (+21 more, 8.8 % together) |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 33225660, 2026-11-17);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 196 clusters | 7.4, 4.58, 4.03, 3.3, 3.27, 2.99, 2.8, 2.75, 2.67, 2.64 ... (+186 more, 63.6 % together) |
| 1FKH | 29 clusters | 28.28, 27.12, 6.98, 6.49, 6.22, 5.87, 5.33, 3.02, 1.83, 1.29 ... (+19 more, 7.5 % together) |
