# 4p6w -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29365313, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `REST2-rotbonds/4p6w-water` (REST2-run1) | `4p6w.input.sdf` |
| 4P6W complex | `REST2-rotbonds/4p6w-complex` | `4p6w.input.sdf` |

Torsions t0..t5 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `24 23 22 17`  #X041-O3 X041-C21 X041-C20 X041-C15 | rotatable
- t1: `26 25 23 22`  #X041-Cl2 X041-C22 X041-C21 X041-C20 | rotatable
- t2: `17 22 27 28`  #X041-C15 X041-C20 X041-O4 X041-C23 | rotatable
- t3: `31 30 28 27`  #X041-C25 X041-C24 X041-C23 X041-O4 | rotatable
- t4: `29 28 27 22`  #X041-O5 X041-C23 X041-O4 X041-C20 | rigid ester C(=O)-O
- t5: `9 11 12 44`  #X041-C9 X041-C10 X041-O2 X041-H10 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).
