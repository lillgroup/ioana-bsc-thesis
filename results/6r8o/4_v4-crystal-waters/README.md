# 6r8o -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365318, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/6r8o-water` (REST2-run1) | `6r8o.input.sdf` |
| 6R8O complex | `REST2-rotbonds/6r8o-complex` | `6r8o.input.sdf` |

Torsions t0..t9 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 20 1 2`  #X081-C1 X081-S1 X081-C2 X081-C3 | rotatable
- t1: `9 21 6 1`  #X081-N1 X081-C20 X081-C7 X081-C2 | rotatable
- t2: `9 23 10 25`  #X081-N1 X081-C22 X081-C10 X081-N2 | rotatable
- t3: `11 25 10 23`  #X081-C11 X081-N2 X081-C10 X081-C22 | rotatable
- t4: `11 27 12 28`  #X081-C11 X081-N3 X081-C12 X081-C23 | rotatable
- t5: `13 28 12 27`  #X081-C13 X081-C23 X081-C12 X081-N3 | rotatable
- t6: `10 23 9 8`  #X081-C10 X081-C22 X081-N1 X081-C9 | rigid amide-like C-N
- t7: `10 25 11 26`  #X081-C10 X081-N2 X081-C11 X081-O2 | rigid amide-like C-N
- t8: `12 27 11 25`  #X081-C12 X081-N3 X081-C11 X081-N2 | rigid amide-like C-N
- t9: `15 18 33 66`  #X081-C15 X081-C18 X081-O3 X081-H32 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
