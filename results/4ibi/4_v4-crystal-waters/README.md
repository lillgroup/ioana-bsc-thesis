# 4ibi -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29363277, 2026-11-02, md-tools 0.6.5.2); state 0, stride 5.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/REST2-conformers-v4/4ibi-water` (REST2-run1) | `4ibi_1D8_A.sdf`, `4ibi_1D8_B.sdf` |
| 4IBI complex | `REST2-complex-xtalwater/4ibi` | `4ibi_1D8_A.sdf` |

Torsions t0..t4 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `0 1 3 4`  #L011-O1 L011-C1 L011-C2 L011-C3 | rotatable
- t1: `6 7 8 9`  #L011-C5 L011-C6 L011-N1 L011-C7 | rotatable
- t2: `11 13 14 15`  #L011-C8 L011-C9 L011-C10 L011-O5 | rotatable
- t3: `13 14 16 17`  #L011-C9 L011-C10 L011-C11 L011-C12 | rotatable
- t4: `8 26 27 28`  #L011-N1 L011-C18 L011-C19 L011-C20 | rotatable

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).



## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597268, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 25 clusters | 22.76, 11.13, 10.46, 8.71, 8.19, 7.71, 7.08, 6.31, 5.83, 2.86 ... (+15 more, 9.0 % together) |
| 4IBI | 2 clusters | 99.8, 0.19 |
