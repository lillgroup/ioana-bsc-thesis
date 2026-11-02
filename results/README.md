# Results

`results/<compound>/<version>/` — one folder per compound, one subfolder per simulation version.
The version number is the same for every compound; a compound has only the versions that were run.

| version | what changed | runs (`runs/` / `$THESIS/`) |
|---|---|---|
| `1_original` | md-tools 0.6.2; ligands as drawn (folate neutral COOH), ligand amides not heated; complexes without crystal waters | FOLn-REST2-50ns, FOL-complex-REST2, REST2-PARA, REST2-conformers (+ `-complex`) |
| `2_v4` | REST2 v4 (md-tools b3603e0): Open Babel pH 7.4 protomers (folate −2, 4ibi/4ibk −2, 4j7d +1), ligand amides heated; complexes still without crystal waters | FOL-REST2-v4, PARA-REST2-v4, REST2-conformers-v4 |
| `3_mace` | folate with Shu-Yu's MACE-OMOL bespoke parameters (3H tautomer, −2), REST2 v4 | FOL-REST2-mace |
| `4_v4-crystal-waters` | as `2_v4`, complexes with the chain-A crystal waters kept | REST2-complex-xtalwater |
| `5_v4-crystal-waters-3H` | as `4_…`, folate as the 3H tautomer (H on N3; 3QL3 only) | REST2-complex-xtalwater/FOL-3QL3-3H |

5mwh and 5ov8 keep the same protomer under the v4 pH model, so they have no `2_v4` run: `1_original`
is their version without crystal waters. For the small molecules the crystal-water runs are, so far,
only in the binding-check overview (`overview/binding-check/crystal-waters/`).

Inside a version folder:
- `symfirst/` — Shu-Yu's symmetry-first clustering: `*_clusters_overview_symfirst-nosingle.png`
  (`allow_single_cluster=False`, **the ones to use**) and `*_symfirst.png` (Shu-Yu's default; folate
  collapses to one cluster). `1_original/` also has the earlier `*_clusters_overview.png`,
  `*_group_torsions.png`, `*_convergence.png` (pre-symmetry-first clustering).
- `bioactive/` — RMSD of each ensemble / cluster to the crystal (bioactive) conformation.
- `replicates/` — REST2 replicates from different starting conformers, pooled and compared.
- folate only: `*_rmsd_vs_crystal*` (ligand conformation vs crystal), `*_water_vs_bound*` /
  `*_water_vs_complexes*` (torsions in water vs bound), `fragments_<PDB>.txt` (pocket-frame drift of
  pteridine / pABA / glutamate and the 3QL3 Asp27 / Arg57 / NADP+ contacts, 5 ns blocks).

`overview/` — figures over several compounds: `binding-check/no-crystal-waters/` and
`binding-check/crystal-waters/` (pocket-frame RMSD, ligand-pocket COM, contacts, native Q per run;
`analysis/binding-check/check_bound.py`), and the README of the 2026-10-15 symmetry-first round.

Scripts that collect several compounds at once (`analysis/run_bioactive.sh`,
`analysis/run_symfirst_0651.sh`) write to `results/_incoming/<tag>/`; sort new outputs from there into
`<compound>/<version>/`.

Details of the `1_original` analyses (same settings for all): state 0, stride 5 (paracetamol 1),
5 ns convergence windows (paracetamol 0.5 ns); 4ibi/4ibk/4j7d/5mwh/5ov8 water 60 ns (run1 + ext1),
folate water 100 ns (FOLn-REST2-50ns), paracetamol 3.5 ns (REST2-PARA, interrupted). Several crystal
lines per torsion mean atoms the crystal cannot tell apart (carboxyl/sulfonyl O, ring flips); O–H
torsions have no line (no H in crystal structures).
