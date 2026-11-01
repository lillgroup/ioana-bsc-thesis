# REST2-rotbonds: 15 new ligands with 4, 5, 6 and 7 rotatable bonds (set up 2026-11-01)

Goal: 5 compounds per rotatable-bond count (4, 5, 6, 7), each in water and bound, with the
small-molecule v4 protocol. Together with the earlier five (4j7d, 4ibi, 4ibk: 5; 5mwh, 5ov8: 6),
15 new ones complete the series:

| rot. bonds | earlier | new (this folder) |
|---|---|---|
| 4 | – | 6cjr, 1v2k, 6cn5, 4p6w, 4k9y |
| 5 | 4j7d, 4ibi, 4ibk | 6s56, 5ufr |
| 6 | 5mwh, 5ov8 | 6r8o, 3b67, 4ase |
| 7 | – | 4erf, 6fo5, 1rd4, 4crf, 1eld |

`compounds.csv`: per compound the PDB entry, deposited ligand code, chain, package residue name
(X01-X15), protein, residues, resolution, affinity, crystal waters and why it was chosen.

## Selection (`selection/`)

The original criteria (`../../filter.py`): MW ~500 (450-550), neutral, <= 40 heavy atoms, smallest
proteins first. Changes, decided 2026-11-01:
- **Rotatable bonds on heavy atoms** (RDKit on the heavy-atom graph, the usual definition).
  `filter.py` counted on the SDF with explicit H, which also counts OH/NH rotors; with that count
  4ibi and 4ibk had 6, they have 5 here.
- "Neutral" is the formal charge of the PDBbind SDF, which is almost always written neutral; the
  packages are re-protonated at pH 7.4 anyway (4ibi, 4ibk were simulated at −2, 4j7d at +1).
- No metal or cofactor in the pocket, one protein chain, no sugars or macrocycles, one entry per
  protein where possible (the refined set has many repeats: VP35, trypsin, Eph kinases, ...).
- **No residues missing inside the chain**: md-tools refuses such gaps (it does not model loops).
  This replaced 4nh7 -> 6cjr, 4p5z -> 4k9y, 6fni -> 5ufr, 5fck -> 4ase; 4tpw -> 4erf because eIF4E
  also binds an m7GTP cap analogue that would need its own parameters.

`selection/scan_refined_set.py` (properties of all refined-set ligands ->
`refined_set_properties.csv`), `selection/screen_candidates.py` (filter + structure screen of each
4-7-bond candidate: ligand chain, gaps, hetero groups, packing contacts -> `candidates_screened.csv`).
The refined set itself is local only (not in git).

## Inputs

- `crystal/<PDB>.pdb` from RCSB; `ccd/<code>_ideal.sdf` (wwPDB Chemical Component Dictionary);
  `pdbbind/<name>_ligand.sdf` (PDBbind's ligand files).
- `scripts/make_inputs.py` -> `packages/<name>.input.sdf` + `.parameterize.config`: deposited heavy-atom
  coordinates with bond orders from PDBbind when its chemistry equals the CCD's (standard InChI), else
  from the CCD -- **4p6w**: PDBbind has mometasone furoate's C20 ketone as an alcohol, the CCD is used.
  Neutral; Open Babel sets the pH 7.4 protomer during parameterisation (Sage + AM1-BCC), the v4
  route of REST2-conformers-v4. These are NEW compounds, so new packages (`packages/<name>/`, made
  by build.sbatch); once made they are reused like every other package, never re-parameterised.
- `scripts/prep.py` -> `<name>-complex/prep/complex.pdb`: the ligand's chain, ligand heavy atoms,
  that chain's crystal waters and ions (md-tools keeps both), altloc A, deposited H dropped.
  Left out, each listed in the script with its reason: crystallisation / cryo additives (SO4, PO4,
  GOL, EDO, PG4, ACY -- as approved for REST2-complex-xtalwater), and one water in 1ELD (HOH 651,
  2.06 A from a ligand carbon: a deposition error). Kept: the structural Ca2+ of trypsin (1v2k,
  23 A from the ligand) and elastase (1eld, 21 A), Cl- in 4crf (13 A).
  Notes: 1rd4 has no deposited waters, 5ufr one; 4p6w's co-activator peptide (chain B) does not
  touch the ligand and is not kept (chain A only, as for every other entry).
- `scripts/make_configs.py` -> per compound
  `<name>-water/` = REST2-conformers-v4/4ibi-water (6 replicas, tau 0-0.4, 60 ns per state), and
  `<name>-complex/` = REST2-complex-xtalwater/4ibi (ff14SB, PROPKA pH 7.4, TIP3P, 1.2 nm cube,
  ligand-only REST2, 6 replicas, tau 0-0.4, 40 ns). The complex mask `:LIGRES` is filled in by
  build.sbatch with the built residue number.

Everything above was run locally on 2026-11-01; its outputs are committed. Nothing was built or
simulated yet: there is no openmm-env-protonation on the local machine.

## On sciCORE

```bash
cd $THESIS/bsc-thesis && git pull
cp -a runs/REST2-rotbonds $THESIS/REST2-rotbonds && cd $THESIS/REST2-rotbonds
bash scripts/submit.sh            # build (array 0-14) -> tests (300 exchanges) -> productions
```

`scripts/submit.sh` chains, per compound: `build.sbatch` task -> water and complex test runs
(`rest2.sbatch`, a100-6hours) -> the productions (a100-1day, auto-resubmitting until done), each
step only if the previous one succeeded. Job ids in `logs/jobs.txt`. Subsets:
`bash scripts/submit.sh 6cjr 4erf`; already built: `NO_BUILD=1 bash scripts/submit.sh 6cjr`;
other GPUs: `GPU=h200 bash scripts/submit.sh`. To build only (and look before running):
`sbatch scripts/build.sbatch`, then `NO_BUILD=1 bash scripts/submit.sh`.

Check after the builds (`logs/build-*.out`): the package's net charge and SMILES after Open Babel
(e.g. 1v2k's aminopyridine and 4erf's COOH should come out charged), that the complex ligand mapped
to the same package, and the mask. After a complex production, rest2.sbatch prints check_run.py's
exchange acceptance and binding report into its log.

## Changes on sciCORE, 2026-11-01

- **6cn5 and 5ufr dropped**: md-tools refused their packages (its parameterisation self-check found the
  package's torsions differ from the openmmforcefields reference: 6cn5 two proper torsions missing,
  5ufr improper torsions in the other cyclic atom order). Replaced by **1l8g** (4 bonds, PTP1B, X16,
  net -2 after Open Babel) and **6qlp** (5 bonds, galectin-3, X17, neutral; PEG fragment left out,
  Cl- kept), from `selection/candidates_screened.csv`. Their inputs come from the CCD (no PDBbind
  files on sciCORE; `make_inputs.py` falls back to the CCD) via `scripts/new_inputs_20261101.sbatch`.
  `make_inputs.py` / `make_configs.py` now take compound names, so earlier rows are not rewritten.
- **1v2k as the pyridinium cation** (user's decision): Open Babel's pH model leaves its
  4-piperidinylpyridine (pKa ~9.5) neutral. `packages/1v2k.input.sdf` is now the cation
  (`1v2k.input.neutral.sdf` the old input). md-tools refuses the cation's package with the same
  improper-order check as 5ufr, also for three other atom orders of the same molecule
  (`scripts/trial_1v2k_atom_order.sbatch`). **1v2k is not running**; it needs a decision (neutral
  package, which passed, or wait for an md-tools fix). Its neutral builds are kept in
  `$THESIS/REST2-rotbonds/_superseded/` (never started).
- Running: 14 compounds, water and complex (`logs/jobs.txt`). The A100 queue was full, so all tests
  and the productions of 6cjr, 4p6w, 4k9y, 6s56, 6r8o, 3b67, 4ase were moved to H200 (h200-6hours,
  2 jobs at a time); the other 7 compounds run on A100.
