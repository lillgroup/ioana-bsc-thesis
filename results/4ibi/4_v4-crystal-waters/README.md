# 4ibi -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620718, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

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

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597268, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 25 clusters | 22.76, 11.13, 10.46, 8.71, 8.19, 7.71, 7.08, 6.31, 5.83, 2.86 ... (+15 more, 9.0 % together) |
| 4IBI | 2 clusters | 99.8, 0.19 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857805, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 17 clusters | 22.76, 11.56, 10.46, 8.71, 8.19, 7.71, 7.08, 6.31, 6.27, 2.97 ... (+7 more, 8.0 % together) |
| 4IBI | 2 clusters | 99.8, 0.19 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857811, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 15 clusters | 22.76, 12.38, 11.56, 8.71, 8.19, 7.71, 7.08, 6.31, 6.27, 2.97 ... (+5 more, 6.1 % together) |
| 4IBI | 2 clusters | 99.8, 0.19 |
