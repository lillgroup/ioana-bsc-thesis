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
