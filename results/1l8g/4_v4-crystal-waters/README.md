# 1l8g -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365392, 2026-11-03, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/1l8g-water` (REST2-run1) | `1l8g.input.sdf` |
| 1L8G complex | `REST2-rotbonds/1l8g-complex` | `1l8g.input.sdf` |

Torsions t0..t6 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `10 9 8 3`  #X161-O2 X161-C8 X161-C7 X161-C3 | rotatable
- t1: `13 12 7 6`  #X161-C9 X161-N1 X161-C6 X161-S1 | rotatable
- t2: `15 14 13 12`  #X161-O4 X161-C10 X161-C9 X161-N1 | rotatable
- t3: `19 18 1 0`  #X161-O7 X161-C11 X161-C1 X161-O1 | rotatable
- t4: `20 19 18 1`  #X161-C12 X161-O7 X161-C11 X161-C1 | rotatable
- t5: `21 20 19 18`  #X161-C13 X161-C12 X161-O7 X161-C11 | rotatable
- t6: `14 13 12 7`  #X161-C10 X161-C9 X161-N1 X161-C6 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
