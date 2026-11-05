# 4ase -- 4_v4-crystal-waters

Made by `analysis/report/report.sbatch` (Slurm job 29620731, 2026-11-04, md-tools 0.6.5.2); state 0, stride 5.

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

Figures: `*_water_clusters.png` / `*_<PDB>_clusters.png` (symmetry-first clustering, no single cluster, recursive splitting while internal RMSD > 2 A, mass floor 0.5 %, 19-of-20 vote), `*_intercluster_rmsd.png` and `*_cluster_network.png`, `*_rmsd_density.png` (ligand conformation vs crystal, water vs complex), `*_binding.png` (complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q), `*_water_vs_complex.png` (torsion distributions).

## Torsion partitioning (barrier cuts)

`*_clusters_cut.png` (+ `_cut_split_tree.png`, `_cut_free_energy.png`, `_cut_groups.csv`): Shu-Yu's
`cut_clustering` (MD-tools main 306b109), `analysis/report/report_cut.sbatch` (Slurm job 29597284, 2026-11-04);
MD-tools defaults (kT_cut 1, mass_floor 0.001, min_p2p pi/3), every frame, same torsions t0.. as above.

| system | clusters | masses (%) |
|---|---|---|
| water | 27 clusters | 29.91, 23.92, 22.14, 20.08, 0.62, 0.25, 0.24, 0.21, 0.2, 0.19 ... (+17 more, 2.2 % together) |
| 4ASE | 5 clusters | 98.12, 0.97, 0.55, 0.22, 0.14 |

`*_cut_intercluster_rmsd.png`, `*_cut_cluster_network.png`: symmetry-corrected heavy-atom RMSD between /
within the barrier-cut clusters >= 1 % of the frames (smaller ones left out, see the title),
`analysis/report/report_cut_rmsd.sbatch` (`shuyu_clustering/cluster_cut_rmsd.py`).
