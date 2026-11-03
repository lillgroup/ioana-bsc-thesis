# 4k9y -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365314, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4k9y-water` (REST2-run1) | `4k9y.input.sdf` |
| 4K9Y complex | `REST2-rotbonds/4k9y-complex` | `4k9y.input.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `12 13 7 2`  #X051-C9 X051-C10 X051-N4 X051-C2 | rotatable
- t1: `22 21 15 16`  #X051-C17 X051-C16 X051-N5 X051-C12 | rotatable
- t2: `29 28 16 15`  #X051-C23 X051-N7 X051-C12 X051-N5 | rotatable
- t3: `29 31 10 9`  #X051-C23 X051-N8 X051-C7 X051-C6 | rotatable
- t4: `30 29 28 16`  #X051-O1 X051-C23 X051-N7 X051-C12 | rigid amide-like C-N
- t5: `10 31 29 28`  #X051-C7 X051-N8 X051-C23 X051-N7 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29468035, 2026-11-03);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 74 clusters | 12.09, 5.09, 4.48, 4.17, 3.97, 3.93, 3.31, 3.07, 2.86, 2.75 ... (+64 more, 54.3 % together) |
| 4K9Y | 2 clusters | 99.65, 0.35 |
