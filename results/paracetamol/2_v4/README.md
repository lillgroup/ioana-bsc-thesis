# paracetamol -- 2_v4

Made by `analysis/report/report.sbatch` (Slurm job 29620717, 2026-11-04, md-tools 0.6.5.2); state 0, stride 1.

| system | run | crystal reference |
|---|---|---|
| water | `water-runs/PARA-REST2-v4` (REST2-run1) |  |

Torsions t0..t2 (the same in every figure; rotatable, then rigid, then polar-H; atom indices of the water ligand):

- t0: `1 3 4 5`  #TYL1-C2 TYL1-N1 TYL1-C3 TYL1-C4 | rotatable
- t1: `0 1 3 4`  #TYL1-C1 TYL1-C2 TYL1-N1 TYL1-C3 | rigid amide-like C-N
- t2: `6 7 8 17`  #TYL1-C5 TYL1-C6 TYL1-O2 TYL1-H7 | polar-H O-H

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597266, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 11 clusters | 52.52, 46.34, 0.16, 0.16, 0.14, 0.14, 0.12, 0.12, 0.1, 0.1 ... (+1 more, 0.1 % together) |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).

## Torsion partitioning (barrier cuts, kT_cut 2)

`*_clusters_cut_kT2.png` (+ `_cut_kT2_split_tree.png`, `_cut_kT2_free_energy.png`, `_cut_kT2_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857816, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 2 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 7 clusters | 52.68, 46.7, 0.16, 0.14, 0.12, 0.1, 0.1 |

## Torsion partitioning (barrier cuts, kT_cut 3)

`*_clusters_cut_kT3.png` (+ `_cut_kT3_split_tree.png`, `_cut_kT3_free_energy.png`, `_cut_kT3_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 30857817, 2026-11-10);
MD-tools defaults (mass_floor 0.001, min_p2p pi/3) with kT_cut 3 (lowest barrier cut E*; default 1), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 7 clusters | 52.68, 46.7, 0.16, 0.14, 0.12, 0.1, 0.1 |
