# 3b67 -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620730, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/3b67-water` (REST2-run1) | `3b67.input.sdf` |
| 3B67 complex | `REST2-rotbonds/3b67-complex` | `3b67.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `8 7 6 4`  #X091-O1 X091-N1 X091-C4 X091-C2 | rotatable
- t1: `14 13 12 5`  #X091-C8 X091-N2 X091-C7 X091-C3 | rotatable
- t2: `17 16 14 13`  #X091-O4 X091-C9 X091-C8 X091-N2 | rotatable
- t3: `20 19 16 14`  #X091-O5 X091-C11 X091-C9 X091-C8 | rotatable
- t4: `21 20 19 16`  #X091-C12 X091-O5 X091-C11 X091-C9 | rotatable
- t5: `22 21 20 19`  #X091-C13 X091-C12 X091-O5 X091-C11 | rotatable
- t6: `15 14 13 12`  #X091-O3 X091-C8 X091-N2 X091-C7 | rigid amide-like C-N
- t7: `14 16 17 36`  #X091-C8 X091-C9 X091-O4 X091-H5 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597282, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 105 clusters | 8.55, 8.53, 5.92, 5.37, 5.36, 5.02, 4.76, 4.57, 4.43, 3.85 ... (+95 more, 43.7 % together) |
| 3B67 | 2 clusters | 77.73, 22.26 |
