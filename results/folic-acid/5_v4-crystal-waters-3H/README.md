# folic-acid -- 5_v4-crystal-waters-3H

Made by `analysis/report/report.sbatch` (Slurm job 29363278, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/FOL-REST2-v4/water` (REST2-run1) | `fol_4LRH_A.sdf`, `fol_3QL3_A.sdf` |
| 4LRH complex | `REST2-complex-xtalwater/FOL-4LRH` | `fol_4LRH_A.sdf` |
| 3QL3 complex | `REST2-complex-xtalwater/FOL-3QL3-3H` | `fol_3QL3_A.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 1 3 4`  #FOL1-O1 FOL1-C1 FOL1-C2 FOL1-N1 | rotatable
- t1: `1 3 4 5`  #FOL1-C1 FOL1-C2 FOL1-N1 FOL1-C3 | rotatable
- t2: `1 3 27 28`  #FOL1-C1 FOL1-C2 FOL1-C17 FOL1-C18 | rotatable
- t3: `4 5 7 8`  #FOL1-N1 FOL1-C3 FOL1-C4 FOL1-C5 | rotatable
- t4: `9 10 13 14`  #FOL1-C6 FOL1-C7 FOL1-N2 FOL1-C10 | rotatable
- t5: `10 13 14 15`  #FOL1-C7 FOL1-N2 FOL1-C10 FOL1-C11 | rotatable
- t6: `13 14 15 16`  #FOL1-N2 FOL1-C10 FOL1-C11 FOL1-N3 | rotatable
- t7: `3 27 28 29`  #FOL1-C2 FOL1-C17 FOL1-C18 FOL1-C19 | rotatable
- t8: `27 28 29 30`  #FOL1-C17 FOL1-C18 FOL1-C19 FOL1-O5 | rotatable
- t9: `3 4 5 6`  #FOL1-C2 FOL1-N1 FOL1-C3 FOL1-O3 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).


## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29369459, 2026-11-02);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 308 clusters | 8.14, 6.04, 5.76, 5.53, 4.82, 4.61, 4.47, 4.0, 2.27, 1.95 ... (+298 more, 52.4 % together) |
| 4LRH | 30 clusters | 44.18, 43.59, 6.98, 0.9, 0.24, 0.22, 0.22, 0.21, 0.21, 0.2 ... (+20 more, 3.0 % together) |
| 3QL3 | 101 clusters | 20.02, 17.59, 11.15, 9.51, 8.9, 8.2, 4.5, 1.6, 1.56, 1.3 ... (+91 more, 15.6 % together) |
