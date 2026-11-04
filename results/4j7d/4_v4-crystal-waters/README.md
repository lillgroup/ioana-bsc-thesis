# 4j7d -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620720, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/REST2-conformers-v4/4j7d-water` (REST2-run1) | `4j7d_I31_A.sdf` |
| 4J7D complex | `REST2-complex-xtalwater/4j7d` | `4j7d_I31_A.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `1 2 3 4`  #L031-C2 L031-C3 L031-C4 L031-C5 | rotatable
- t1: `5 6 7 8`  #L031-N1 L031-C6 L031-C7 L031-C8 | rotatable
- t2: `7 16 17 18`  #L031-C7 L031-C16 L031-O1 L031-C17 | rotatable
- t3: `16 17 18 19`  #L031-C16 L031-O1 L031-C17 L031-C18 | rotatable
- t4: `3 21 22 23`  #L031-C4 L031-C19 L031-C20 L031-C21 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597270, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 8 clusters | 39.32, 38.51, 6.12, 5.78, 5.15, 4.65, 0.24, 0.23 |
| 4J7D | 10 clusters | 71.83, 12.47, 6.86, 5.73, 1.67, 0.53, 0.37, 0.24, 0.19, 0.11 |
