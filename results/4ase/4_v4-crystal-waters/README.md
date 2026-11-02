# 4ase -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365320, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4ase-water` (REST2-run1) | `4ase.input.sdf` |
| 4ASE complex | `REST2-rotbonds/4ase-complex` | `4ase.input.sdf` |

Torsions t0..t7 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `11 10 9 3`  #X101-C10 X101-O1 X101-C9 X101-C4 | rotatable
- t1: `12 11 10 9`  #X101-C11 X101-C10 X101-O1 X101-C9 | rotatable
- t2: `21 18 13 14`  #X101-C16 X101-N2 X101-C12 X101-C13 | rotatable
- t3: `30 19 5 0`  #X101-C21 X101-O3 X101-C6 X101-C1 | rotatable
- t4: `31 20 0 1`  #X101-C22 X101-O4 X101-C1 X101-C2 | rotatable
- t5: `27 24 23 21`  #X101-C20 X101-C17 X101-N3 X101-C16 | rotatable
- t6: `22 21 18 13`  #X101-O5 X101-C16 X101-N2 X101-C12 | rigid amide-like C-N
- t7: `24 23 21 18`  #X101-C17 X101-N3 X101-C16 X101-N2 | rigid amide-like C-N

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
