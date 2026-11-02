# paracetamol -- 2_v4

Made by `analysis/report/report.sbatch` (Slurm job 29365307, 2026-11-02, md-tools 0.6.5.2); state 0, stride 1.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/PARA-REST2-v4` (REST2-run1) |  |

Torsions t0..t2 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `1 3 4 5`  #TYL1-C2 TYL1-N1 TYL1-C3 TYL1-C4 | rotatable
- t1: `0 1 3 4`  #TYL1-C1 TYL1-C2 TYL1-N1 TYL1-C3 | rigid amide-like C-N
- t2: `6 7 8 17`  #TYL1-C5 TYL1-C6 TYL1-O2 TYL1-H7 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster), `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).



## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29369458, 2026-11-02);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 11 clusters | 52.52, 46.34, 0.16, 0.16, 0.14, 0.14, 0.12, 0.12, 0.1, 0.1 ... (+1 more, 0.1 % together) |
