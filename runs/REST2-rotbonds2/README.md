# REST2-rotbonds2: 28 more ligands, 4-9 rotatable bonds (set up 2026-11-12)

Goal: 8 compounds per heavy-atom rotatable-bond count 4-9, together with REST2-rotbonds (4-7), the
original five (4ibi, 4ibk, 4j7d, 5mwh, 5ov8) and folate (9). Per compound a water run (6 replicas,
tau 0-0.4, 60 ns per state) and a complex run (crystal waters, ff14SB, PROPKA pH 7.4, TIP3P, 1.2 nm cube,
ligand-only REST2, 6 replicas, 40 ns), the settings of REST2-rotbonds; new packages (Open Babel pH 7.4,
Sage + AM1-BCC, ligand amides scaled: REST2 v4).

**Difference to REST2-rotbonds: MD-tools 42a949b (feature/reduce-protonation, env openmm-env-reduce):
Reduce decides the Asn/Gln/His flips and His tautomers of the protein** (`protonation.flips: reduce`).

Selection (`ligand-selection/rotbonds/`): the REST2-rotbonds criteria (neutral, <= 40 heavy atoms, one
chain, no metal/cofactor in the pocket, no B/Si/Se/P, structure screen: no internal chain gaps, no other
hetero groups than additives/ions), MW 450-550; for 4, 8 and 9 rotatable bonds the window was widened to
400-600 (too few usable compounds otherwise; noted per row in compounds.csv). 4-7: candidates already
screened (`candidates_screened.csv`); 4 (wider MW), 8, 9: `screen_by_code.py` (the ligand found by its
residue code, as PDBbind's files are not on sciCORE) -> `candidates_screened_4_8_9.csv`, `_8b`, `_9b`.

Replaced during set-up (each noted in compounds.csv): covalent ligands 1k1i, 1ppc, 4jz1 (ligand atoms
bonded to the enzyme missing); ligands with undeposited atoms 5mks, 4h3f; packages md-tools refused
(parameterisation self-check, as 5ufr / 6cn5 before) 6g2r, 6el5, 6p9e.
prep.py leaves out, with the reason printed: crystallisation additives (now also SCN, PEG) and single
deposited waters overlapping the solute (4jfk 532, 1fkh 124/223, 6eog 1104/1113).

On sciCORE: `$THESIS/REST2-rotbonds2` (scripts/make_all_inputs.sbatch, then `bash scripts/submit.sh`;
job ids in logs/jobs.txt). This folder holds the inputs only.
