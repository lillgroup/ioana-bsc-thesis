# 1eld -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620736, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/1eld-water` (REST2-run1) | `1eld.input.sdf` |
| 1ELD complex | `REST2-rotbonds/1eld-complex` | `1eld.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `8 7 6 0`  #X151-C4 X151-C3 X151-N1 X151-C1 | rotatable
- t1: `9 8 7 6`  #X151-O2 X151-C4 X151-C3 X151-N1 | rotatable
- t2: `11 10 7 6`  #X151-C6 X151-C5 X151-C3 X151-N1 | rotatable
- t3: `12 11 10 7`  #X151-C7 X151-C6 X151-C5 X151-C3 | rotatable
- t4: `19 18 17 8`  #X151-C13 X151-C12 X151-N2 X151-C4 | rotatable
- t5: `20 19 18 17`  #X151-O3 X151-C13 X151-C12 X151-N2 | rotatable
- t6: `24 23 22 19`  #X151-C16 X151-C15 X151-N3 X151-C13 | rotatable
- t7: `7 6 0 1`  #X151-C3 X151-N1 X151-C1 X151-C2 | rigid amide-like C-N
- t8: `18 17 8 7`  #X151-C12 X151-N2 X151-C4 X151-C3 | rigid amide-like C-N
- t9: `23 22 19 18`  #X151-C15 X151-N3 X151-C13 X151-C12 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597292, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 44 clusters | 38.28, 29.49, 11.48, 6.52, 4.11, 2.25, 0.87, 0.6, 0.56, 0.47 ... (+34 more, 5.4 % together) |
| 1ELD | 2 clusters | 97.09, 2.9 |
