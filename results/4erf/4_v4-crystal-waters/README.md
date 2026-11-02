# 4erf -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365321, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4erf-water` (REST2-run1) | `4erf.input.sdf` |
| 4ERF complex | `REST2-rotbonds/4erf-complex` | `4erf.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `8 7 6 0`  #X111-C8 X111-C7 X111-C6 X111-C1 | rotatable
- t1: `8 9 10 11`  #X111-C8 X111-C9 X111-C10 X111-C11 | rotatable
- t2: `12 11 10 9`  #X111-O1 X111-C11 X111-C10 X111-C9 | rotatable
- t3: `19 18 17 15`  #X111-C15 X111-C14 X111-N1 X111-C13 | rotatable
- t4: `20 19 18 17`  #X111-C16 X111-C15 X111-C14 X111-N1 | rotatable
- t5: `22 21 18 17`  #X111-O4 X111-C17 X111-C14 X111-N1 | rotatable
- t6: `26 25 24 7`  #X111-C21 X111-C20 X111-C19 X111-C7 | rotatable
- t7: `18 21 22 51`  #X111-C14 X111-C17 X111-O4 X111-H20 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
