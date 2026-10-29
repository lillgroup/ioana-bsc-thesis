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

## Tautomer test: FOL-3QL3-3H (2026-10-29)

FOL-3QL3 with folate as the 3H tautomer (H on pteridine N3 instead of N1). In 3QL3, Asp27 OD2 is
2.64 A from N3 (deprotonated Asp27, so N3 must carry the H); the v4 package has N1-H. New package
`packages/FOL-3H` (param_a633e316c1e8; same route as v4: neutral SMILES `FOL-3H.smi`, Open Babel pH 7.4,
net -2, Sage 2.2.1 + AM1-BCC), made by `scripts/param_fol_3h.sbatch`, which checks the 3H dianion
after Open Babel. Built with `build.sbatch --array=7`; everything else as FOL-3QL3.
Chain: param 29052311 -> build 29052312 -> test 29052313 -> production 29052314.
