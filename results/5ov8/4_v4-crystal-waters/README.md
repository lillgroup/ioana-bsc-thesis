# 5ov8 -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620722, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/REST2-conformers/5ov8` (REST2-run1 + REST2-run1-ext1) | `5ov8_AXN_A.sdf`, `5ov8_AXN_B.sdf` |
| 5OV8 complex | `REST2-complex-xtalwater/5ov8` | `5ov8_AXN_A.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `1 2 3 4`  #L051-C2 L051-C3 L051-S1 L051-O1 | rotatable
- t1: `2 3 6 7`  #L051-C3 L051-S1 L051-N1 L051-C4 | rotatable
- t2: `3 6 7 8`  #L051-S1 L051-N1 L051-C4 L051-C5 | rotatable
- t3: `7 8 9 10`  #L051-C4 L051-C5 L051-N2 L051-C6 | rotatable
- t4: `0 28 29 30`  #L051-C1 L051-C20 L051-C21 L051-C22 | rotatable
- t5: `28 29 30 31`  #L051-C20 L051-C21 L051-C22 L051-C23 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29598920, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 72 clusters | 7.62, 7.58, 7.07, 6.41, 6.27, 6.24, 5.19, 5.12, 3.12, 3.05 ... (+62 more, 42.3 % together) |
| 5OV8 | 27 clusters | 20.27, 18.46, 18.16, 16.67, 8.58, 4.76, 3.15, 2.48, 2.23, 2.23 ... (+17 more, 3.0 % together) |
