# Folate RMSD to the crystal conformation: complexes with crystal waters

From `runs/REST2-complex-xtalwater/scripts/rmsd_vs_crystal.sbatch` (sciCORE job 29045718,
2026-10-29). Folate superposed on folate (shape, not position in the pocket; for that see
`results/binding-check/xtalwater/`), REST2 state 0, best over symmetry-equivalent atoms. Folic acid
in water (reference) is read from `results/symfirst-2026-10-15/1_v4-new-rest2/folic-acid-v4_rmsd_vs_crystal.csv`.

| run | median RMSD | < 1 A | < 1.5 A |
|---|---|---|---|
| 4LRH, crystal waters | 0.63 A | 90 % | 97 % |
| 3QL3, crystal waters | 1.98 A | 17 % | 23 % |
| in water vs 4LRH crystal | 2.13 A | 0 % | 7 % |
| in water vs 3QL3 crystal | 2.54 A | 0 % | 4 % |

The `.csv` holds the new runs only (columns run, crystal, time_ns, rmsd_A, as the earlier file).
