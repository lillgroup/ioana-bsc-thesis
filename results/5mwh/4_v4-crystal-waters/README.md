# 5mwh -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620721, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/REST2-conformers/5mwh` (REST2-run1 + REST2-run1-ext1) | `5mwh_UWX_A.sdf`, `5mwh_UWX_B.sdf` |
| 5MWH complex | `REST2-complex-xtalwater/5mwh` | `5mwh_UWX_A.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `1 2 3 4`  #L041-C2 L041-C3 L041-S1 L041-O1 | rotatable
- t1: `2 3 6 7`  #L041-C3 L041-S1 L041-N1 L041-C4 | rotatable
- t2: `3 6 7 8`  #L041-S1 L041-N1 L041-C4 L041-C5 | rotatable
- t3: `7 8 9 10`  #L041-C4 L041-C5 L041-N2 L041-C6 | rotatable
- t4: `0 29 30 31`  #L041-C1 L041-C20 L041-C21 L041-C22 | rotatable
- t5: `29 30 31 32`  #L041-C20 L041-C21 L041-C22 L041-C23 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29598919, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 50 clusters | 9.18, 8.76, 7.8, 7.75, 7.71, 6.44, 6.14, 5.76, 4.09, 3.97 ... (+40 more, 32.4 % together) |
| 5MWH | 19 clusters | 31.4, 15.32, 10.6, 10.19, 9.03, 8.01, 6.52, 6.27, 1.02, 0.23 ... (+9 more, 1.4 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
