# 1l8g -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620726, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

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

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597277, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 17 clusters | 20.58, 17.35, 14.37, 10.38, 9.73, 8.42, 7.51, 6.95, 1.9, 1.53 ... (+7 more, 1.3 % together) |
| 1L8G | 7 clusters | 82.85, 8.39, 7.64, 0.64, 0.18, 0.15, 0.15 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857824, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 11 clusters | 31.72, 20.58, 16.69, 15.93, 10.38, 3.43, 0.33, 0.24, 0.24, 0.24 ... (+1 more, 0.2 % together) |
| 1L8G | 2 clusters | 99.85, 0.15 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857825, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 9 clusters | 42.1, 20.58, 20.12, 15.93, 0.33, 0.24, 0.24, 0.24, 0.22 |
| 1L8G | 2 clusters | 99.85, 0.15 |
