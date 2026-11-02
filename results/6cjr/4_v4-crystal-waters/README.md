# 6cjr -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365312, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/6cjr-water` (REST2-run1) | `6cjr.input.sdf` |
| 6CJR complex | `REST2-rotbonds/6cjr-complex` | `6cjr.input.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 6 1 2`  #X011-C1 X011-N2 X011-C2 X011-C3 | rotatable
- t1: `1 6 0 14`  #X011-C2 X011-N2 X011-C1 X011-C10 | rotatable
- t2: `4 8 7 10`  #X011-C4 X011-C5 X011-N3 X011-N4 | rotatable
- t3: `3 26 2 1`  #X011-N1 X011-C17 X011-C3 X011-C2 | rotatable
- t4: `11 12 25 52`  #X011-C7 X011-C8 X011-O3 X011-H20 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
