# 1rd4 -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365406, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/1rd4-water` (REST2-run1) | `1rd4.input.sdf` |
| 1RD4 complex | `REST2-rotbonds/1rd4-complex` | `1rd4.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `3 2 1 0`  #X131-C2 X131-C1 X131-N1 X131-O1 | rotatable
- t1: `9 8 6 5`  #X131-C8 X131-C7 X131-C5 X131-C4 | rotatable
- t2: `14 13 12 11`  #X131-C11 X131-N3 X131-C10 X131-N2 | rotatable
- t3: `24 23 3 2`  #X131-C18 X131-S1 X131-C2 X131-C1 | rotatable
- t4: `25 24 23 3`  #X131-C19 X131-C18 X131-S1 X131-C2 | rotatable
- t5: `31 30 25 24`  #X131-C24 X131-O3 X131-C19 X131-C18 | rotatable
- t6: `32 31 30 25`  #X131-C25 X131-C24 X131-O3 X131-C19 | rotatable
- t7: `20 19 16 15`  #X131-O2 X131-C15 X131-N4 X131-C12 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
