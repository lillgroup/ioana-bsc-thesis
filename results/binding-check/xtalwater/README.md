# Binding check: complex REST2 with crystal waters

`check_bound.py` (one folder up) on the 8 runs in `runs/REST2-complex-xtalwater` (v4 REST2, chain-A
crystal waters kept), state 0. Run on sciCORE with `runs/REST2-complex-xtalwater/scripts/analyse_3h.sbatch`
(2026-10-31, Slurm job 29309433; the 7 first runs unchanged from job 29028590); the older runs are not on that cluster, so they are not in these figures.

| run | RMSD pocket frame, last 10 % | Q, last 10 % | verdict |
|---|---|---|---|
| FOL-4LRH | 1.6 A | 0.76 | bound, crystal pose |
| FOL-3QL3 | 4.9 A | 0.15 | bound, crystal pose lost (switch at ~10-11 ns) |
| FOL-3QL3-3H (folate 3H tautomer) | 2.7 A | 0.72 | bound, crystal pose (pteridine 0.5-0.7 A, Asp27 2.8 A throughout; glutamate tail leaves Arg57) |
| 4ibi | 1.0 A | 0.76 | bound, crystal pose |
| 4ibk | 1.1 A | 0.83 | bound, crystal pose |
| 4j7d | 1.5 A | 0.62 | bound, crystal pose |
| 5mwh | 3.3 A | 0.60 | bound, rearranged |
| 5ov8 | 2.3 A | 0.63 | bound, crystal pose |

No ligand dissociated (closest ligand-receptor contact never above 3.6 A).
