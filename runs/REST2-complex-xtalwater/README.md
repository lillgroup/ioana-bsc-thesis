# REST2 complexes with crystal waters (2026-10-26)

Reruns of the 7 complex REST2 runs: FOL-4LRH, FOL-3QL3 (from FOL-REST2-v4), 4ibi, 4ibk, 4j7d
(from REST2-conformers-v4) and 5mwh, 5ov8 (from REST2-conformers-complex). The old runs started
from prep.py output without crystal waters; 3QL3 folate drifted out of the pocket. Here the
chain-A crystal waters are kept, as md-tools does for a deposited structure. Packages, configs and
protein/ligand coordinates are unchanged.

Inputs not included, approved by Shu-Yu on 2026-10-26:

| PDB | not included | reason |
|---|---|---|
| 4IBK | 2 DMS (one 3.2 A from the ligand), 1 GOL | crystallisation solution |
| 4J7D | 2 SO4 | crystallisation solution |
| 5OV8 | 1 NO3 | crystallisation solution |
| 5MWH | 2 NO3 (chain B) | crystallisation solution |
| 3QL3 | 2 Mn2+ (15-25 A from FOL/NAP; Asp70, Asp116/His149) | would need parameters; far from the site |
| 4LRH | NAG on Asn47, Asn139, Asn179 (glycosylated Asn -> Asn; nearest 8.8 A from FOL) | would need parameters; far from the site |

Also: chain A only (biological assembly 1 is the chain-A monomer for all seven entries), altloc A.

Scripts: `scripts/prep.py` (crystal -> prep/complex.pdb with waters; every left-out group is listed
there with its reason), `scripts/build.sbatch` (array 0-6), `scripts/rest2.sbatch`,
`scripts/check_run.py`. Test jobs (300 exchanges) are listed in `logs/test-jobs.txt`.
