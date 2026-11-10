# 4p6w -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620724, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4p6w-water` (REST2-run1) | `4p6w.input.sdf` |
| 4P6W complex | `REST2-rotbonds/4p6w-complex` | `4p6w.input.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `24 23 22 17`  #X041-O3 X041-C21 X041-C20 X041-C15 | rotatable
- t1: `26 25 23 22`  #X041-Cl2 X041-C22 X041-C21 X041-C20 | rotatable
- t2: `17 22 27 28`  #X041-C15 X041-C20 X041-O4 X041-C23 | rotatable
- t3: `31 30 28 27`  #X041-C25 X041-C24 X041-C23 X041-O4 | rotatable
- t4: `29 28 27 22`  #X041-O5 X041-C23 X041-O4 X041-C20 | rigid ester C(=O)-O
- t5: `9 11 12 44`  #X041-C9 X041-C10 X041-O2 X041-H10 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597275, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 11 clusters | 56.62, 27.95, 4.08, 3.42, 3.05, 2.29, 1.61, 0.44, 0.22, 0.17 ... (+1 more, 0.1 % together) |
| 4P6W | 5 clusters | 80.23, 9.87, 8.44, 1.04, 0.43 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857820, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 8 clusters | 87.99, 6.37, 3.05, 1.61, 0.44, 0.22, 0.17, 0.15 |
| 4P6W | 3 clusters | 90.09, 9.48, 0.43 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857821, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 6 clusters | 87.99, 6.37, 4.83, 0.44, 0.22, 0.15 |
| 4P6W | 3 clusters | 90.09, 9.48, 0.43 |
