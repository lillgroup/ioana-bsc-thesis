# 4ibk -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365308, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/REST2-conformers-v4/4ibk-water` (REST2-run1) | `4ibk_1DE_A.sdf`, `4ibk_1DE_B.sdf` |
| 4IBK complex | `REST2-complex-xtalwater/4ibk` | `4ibk_1DE_A.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 1 3 4`  #L021-O1 L021-C1 L021-C2 L021-C3 | rotatable
- t1: `6 7 8 9`  #L021-C5 L021-C6 L021-N1 L021-C7 | rotatable
- t2: `11 13 14 15`  #L021-C8 L021-C9 L021-C10 L021-O5 | rotatable
- t3: `13 14 16 17`  #L021-C9 L021-C10 L021-C11 L021-S1 | rotatable
- t4: `8 22 23 24`  #L021-N1 L021-C15 L021-C16 L021-C17 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).



## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597269, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 26 clusters | 17.41, 15.26, 10.24, 9.13, 8.89, 6.82, 6.42, 6.18, 4.72, 3.37 ... (+16 more, 11.6 % together) |
| 4IBK | 2 clusters | 93.67, 6.33 |
