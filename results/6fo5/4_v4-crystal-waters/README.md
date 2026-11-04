# 6fo5 -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365390, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/6fo5-water` (REST2-run1) | `6fo5.input.sdf` |
| 6FO5 complex | `REST2-rotbonds/6fo5-complex` | `6fo5.input.sdf` |

Torsions t0..t6 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `10 9 4 12`  #X121-C11 X121-C10 X121-C5 X121-N2 | rotatable
- t1: `18 12 4 9`  #X121-S1 X121-N2 X121-C5 X121-C10 | rotatable
- t2: `12 18 7 6`  #X121-N2 X121-S1 X121-C8 X121-C7 | rotatable
- t3: `7 18 12 4`  #X121-C8 X121-S1 X121-N2 X121-C5 | rotatable
- t4: `14 24 21 1`  #X121-N3 X121-C18 X121-C15 X121-C2 | rotatable
- t5: `21 24 14 25`  #X121-C15 X121-C18 X121-N3 X121-C19 | rotatable
- t6: `36 35 32 27`  #X121-C25 X121-C24 X121-N6 X121-C21 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).


## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597288, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 180 clusters | 6.97, 5.9, 4.31, 3.5, 3.27, 2.91, 2.8, 2.47, 2.38, 2.31 ... (+170 more, 63.2 % together) |
| 6FO5 | 52 clusters | 46.36, 9.75, 9.63, 9.05, 5.05, 4.96, 1.77, 1.75, 1.18, 1.06 ... (+42 more, 9.4 % together) |
