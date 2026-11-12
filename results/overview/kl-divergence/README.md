# KL and Jensen-Shannon divergence: bound vs free ligand ensemble

Barrier-cut (torsion partitioning) clusters, REST2 state 0. Scripts and formulas: `analysis/report/README.md`.

| figure | what |
|---|---|
| `kl_cross_cut.png` | cross KL (Shu-Yu): y = KL(pi^p \|\| pi^w) on the water clusters (complex frames assigned), x = KL(pi^w \|\| pi^p) on the complex clusters (water frames assigned); colour = rotatable bonds |
| `kl_cross_cut_regions.png` | the same with the lock and key / conformational selection / induced fit regions |
| `kl_cross_cut_examples.png`, `_examples_regions.png` | only the four corner cases 4j7d, 4crf, 4ase, 1l8g |
| `kl_joint_cut.png` | KL in both directions on ONE joint clustering of water + complex (x axis dominated by the 0.5-frame pseudocount, see `*`) |
| `js_joint_halves_cut.png` | Jensen-Shannon halves KL(pi^w \|\| m), KL(pi^p \|\| m), m = (pi^w + pi^p)/2, each in [0, ln 2] |
| `js_joint_halves_cut_regions.png` | the same with the lock and key / conformational selection / induced fit regions: the cross-KL rules in JS-half units, y = ln(2/(1+w)) for a bound conformer holding a fraction w of water (right axis) |
| `js_joint_ranked_cut.png` | JS per ligand, ranked |

- `data/`: every number behind the figures (`kl_cross_cut.csv`, `kl_joint_cut.csv` with all compounds; one CSV per compound and method).
- `other/`: HDBSCAN versions (noise frames dropped, not used for interpretation) and the superseded pooled KL (`kl_divergence_*.png`).
