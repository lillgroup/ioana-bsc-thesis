# Binding check: complex REST2 without crystal waters

`check_bound.py` (one folder up) on the 14 complex runs built from `prep/complex.pdb` files that
DROPPED every crystal water (`FOL-complex-REST2/scripts/prep.py`, `DROP = {"HOH", ...}`; the
small-molecule preps likewise), state 0, every 10th frame. Run locally (WSL) on 2026-10-22 on the
trajectories copied from sbrinz. Compare with `../xtalwater/` (the same check on the reruns with
the chain-A crystal waters kept).

Measures (see `check_bound.py`): ligand heavy-atom RMSD to the crystal pose in the pocket frame,
ligand-pocket centre-of-mass distance, all ligand-receptor contacts (< 4 A; receptor = protein +
cofactor, NADP+ in 3QL3) and native contacts Q (crystal pairs < 4 A still < 4.5 A).

| run | RMSD pocket frame, last 10 % | Q, last 10 % | verdict |
|---|---|---|---|
| FOL-REST2-v4/4LRH | 1.2 A | 0.83 | bound, crystal pose |
| FOL-REST2-v4/3QL3 | 7.8 A | 0.12 | bound, crystal pose lost |
| FOL-REST2-mace/4LRH | 3.1 A | 0.33 | bound, rearranged |
| FOL-REST2-mace/3QL3 | 5.7 A | 0.65 | bound, rearranged |
| FOL-complex-REST2/4LRH | 1.6 A | 0.73 | bound, crystal pose |
| FOL-complex-REST2/3QL3 | 8.3 A | 0.23 | bound, crystal pose lost |
| REST2-conformers-v4/4ibi-complex | 1.0 A | 0.77 | bound, crystal pose |
| REST2-conformers-v4/4ibk-complex | 1.1 A | 0.85 | bound, crystal pose |
| REST2-conformers-v4/4j7d-complex | 1.4 A | 0.57 | bound, rearranged |
| REST2-conformers-complex/4ibi | 1.1 A | 0.76 | bound, crystal pose |
| REST2-conformers-complex/4ibk | 1.8 A | 0.71 | bound, crystal pose |
| REST2-conformers-complex/4j7d | 2.3 A | 0.41 | bound, rearranged |
| REST2-conformers-complex/5mwh | 3.3 A | 0.61 | bound, rearranged |
| REST2-conformers-complex/5ov8 | 2.5 A | 0.60 | bound, crystal pose |

No ligand dissociated (closest ligand-receptor contact never above 3.6 A). In FOL-REST2-v4/3QL3
the state-0 trajectory mixes replicas: one (walker 6) has the folate out of the pocket
(ligand-pocket COM ~10.5 A vs 2.1 A in the crystal), the others keep it bound.

`com_filter/` (`com_filter.py`, one folder up): per-frame ligand-pocket COM distance for every run,
bound if <= 6 A. `com_histograms.png` and `com_fraction_bound.csv` are in git; the per-frame
`*_com.csv` files (frame, time, walker, COM, bound; 7 MB) are kept locally only.
