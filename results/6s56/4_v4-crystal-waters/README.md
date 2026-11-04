# 6s56 -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365316, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/6s56-water` (REST2-run1) | `6s56.input.sdf` |
| 6S56 complex | `REST2-rotbonds/6s56-complex` | `6s56.input.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 5 3 2`  #X061-C1 X061-S1 X061-N1 X061-C3 | rotatable
- t1: `3 5 0 8`  #X061-N1 X061-S1 X061-C1 X061-C5 | rotatable
- t2: `11 10 9 8`  #X061-C7 X061-N2 X061-C6 X061-C5 | rotatable
- t3: `10 11 1 13`  #X061-N2 X061-C7 X061-C2 X061-N3 | rotatable
- t4: `14 13 1 11`  #X061-C8 X061-N3 X061-C2 X061-C7 | rotatable
- t5: `1 11 10 9`  #X061-C2 X061-C7 X061-N2 X061-C6 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).


## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597278, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 26 clusters | 12.87, 11.92, 11.59, 10.6, 10.1, 9.79, 5.7, 5.46, 5.17, 5.08 ... (+16 more, 11.7 % together) |
| 6S56 | 3 clusters | 99.31, 0.57, 0.11 |
