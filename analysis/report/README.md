# Per-compound report: the figures for the thesis

One Slurm job per compound (`report.sbatch <compound>`, or `bash submit.sh [compound ...]`) makes,
from the compound's current water and complex runs (`compounds.csv`):

| figure (`results/<compound>/<version>/`) | what | script |
|---|---|---|
| `<c>_water_clusters.png`, `<c>_<PDB>_clusters.png` | Shu-Yu's symmetry-first clustering, `--no-single-cluster`, water and each complex | `../shuyu_clustering/cluster_prune.py` |
| `<c>_rmsd_density.png` | ligand conformation RMSD to the crystal, water vs complex (density) | `rmsd_density.py` |
| `<c>_binding.png` | complex: pocket-frame RMSD, ligand-pocket COM, contacts, native Q over time | `../binding-check/check_bound.py --runs` |
| `<c>_water_vs_complex.png` | torsion distributions, water vs complex, one panel per torsion | `water_vs_complex.py` |

plus the CSVs behind them and a `README.md` (runs, crystal references, torsion list, md-tools version,
Slurm job). `prepare.py` writes the inputs to `$THESIS/analysis-report/<compound>/` (outside git).

**One torsion numbering per compound.** `prepare.py` runs `select_torsions.py` once on the water
ligand (rotatable, then rigid, then polar-H; ring torsions are not used) and maps those torsions onto
every complex ligand by matching the heavy-atom graph. t<k> is therefore the same bond in all four
figures, also when the complex uses a different package (3QL3: folate as the 3H tautomer).

`compounds.csv`: water run (+ extensions joined, e.g. 5mwh/5ov8 run1 + ext1), complexes as
`<PDB>=<run>=<crystal SDF>`, crystal SDFs for the water clustering, stride, automorphism cap, and `gauge`:
the torsion position for cluster_prune's --gauge where one symmetry moves two torsions together (first
torsion of that pair; folate t3 = the pABA ring, paracetamol t0 = aryl as in Shu-Yu's tutorial). Run folders
are relative to `$THESIS`, crystal SDFs to the repository. Environment: `analysis-env-main`
(MD-tools c30c1c0 = 0.6.5.2).

A complex whose clustering finds no cluster (one conformation for the whole run, e.g. 4ibi bound)
gets `<c>_<PDB>_NO_CLUSTER.txt` instead of a figure: only `--no-single-cluster` is run.

## Second clustering: torsion partitioning (barrier cuts)

`report_cut.sbatch <compound>` runs `../shuyu_clustering/cluster_cut.py` -- Shu-Yu's `cut_clustering`
(MD-tools main 306b109; docs: openmm_methods/analysis/torsional-clustering/torsion-partitioning) -- on
the same prepared inputs, water and every complex, MD-tools' defaults, every frame. Figures:
`<c>_<system>_clusters_cut.png` (same layout and t-numbering as the symmetry-first figures; clusters
>= 1 % drawn, the rest summed in the title; each cluster in the labelling nearest its own central frame),
plus MD-tools' `_cut_split_tree.png` and `_cut_free_energy.png`, and `_cut_groups.csv`. It appends a
section to the folder's README (rerun it after `report.sbatch`, which rewrites the README).
Environment: `analysis-env-cut` (a clone of `analysis-env-main` with md-tools from the `MD-tools-cut`
worktree at 306b109), so `analysis-env-main` and the symmetry-first results stay on c30c1c0.

## Conformational entropy across ligands

`conformational_entropy.py` (`sbatch conformational_entropy.sbatch`): S_c = -sum pi_c ln pi_c (k_B) for
every ligand in water and complex,
from the populations of each clustering method -> `results/overview/conformational-entropy/`
(`_hdbscan.png`, `_cut.png`, `.csv`). Each figure uses one method with one set of settings; the script
checks every run's own record and stops if they differ. HDBSCAN: populations over the assigned frames
(the unassigned fraction is in the CSV); a run with no cluster is drawn at S_c = 0, hollow. Error bars:
block standard error over 5 time blocks. x: rotatable bonds, RDKit on the heavy-atom graph.

## KL divergence: bound vs free ensemble

`kl_divergence.py` (`sbatch kl_divergence.sbatch`): KL(pi^complex || pi^water) on ONE shared clustering
of water + complex frames (equal numbers from each) per (ligand, complex), for both methods, each in its
own environment -> `results/overview/kl-divergence/` (`other/kl_divergence_{cut,hdbscan}.png`, `data/kl_{cut,hdbscan}.csv`).
Next to KL: S_w (KL of selecting a typical water conformer) and the fraction of the bound ensemble in
clusters with no water frame (induced fit). Pseudocount 0.5 frames per cluster; block-SE error bars.

## KL divergence, Shu-Yu's protocol: clusters of one ensemble, frames of the other assigned

`kl_cross.py` (`bash kl_cross_submit.sh`: one 30-min task per compound and method, then the figures):
water clusters + complex frames assigned -> KL(pi^p || pi^w) (y); complex clusters + water frames
assigned -> KL(pi^w || pi^p) (x) -> `results/overview/kl-divergence/kl_cross_cut*.png` (HDBSCAN: `other/`), CSVs in `data/`.
Barrier cuts: a frame is walked down the reference's cut tree (MD-tools has no assignment for this
method). HDBSCAN: THDBSCAN.predict; its unassigned fraction is in the CSV (very large for water frames
on complex clusters, so the HDBSCAN x values rest on few frames). This replaces the pooled version
above (`kl_divergence.py`, one clustering of water + complex together).

`kl_cross_{cut,hdbscan}_regions.png`: the same scatter with the binding regimes as background
(`regime()` in `kl_cross.py`; the CSV's `regime` column), y = KL(p||w) ~ -ln(population of the bound
conformer in water) when the complex sits mostly in one cluster (right axis):
- lock and key: both KL < 1 k_B (`--lock-key-max`): free and bound ensembles nearly the same;
- induced fit: KL(p||w) >= 3 k_B (`--induced-fit-min`): bound conformer <= ~5 % of the free ensemble;
- conformational selection: the rest: the bound conformer is a populated state in water, binding shifts
  the populations towards it.
Thresholds are our choice, not a standard; KL measures populations, not the binding mechanism's order of
events. The HDBSCAN version drops unassigned frames, so its regimes are unreliable (use the cut one).

## Cost of the bioactive conformation in water (Zivanovic et al., JCTC 2020)

`bioactive_cost.py` (`sbatch bioactive_cost.sbatch`): dG_dist = -ln P_b (P_b = population of the
barrier-cut WATER cluster the crystal falls in, the gold frame of `*_water_clusters_cut.png`) and
dG_strain = -ln(P_b/P_0) (P_0 = most populated water cluster), in k_BT, block SE over 5 time blocks;
cluster-free cross-check: fraction of water frames within 1.0 / 1.5 / 2.0 A symmetry-corrected
heavy-atom RMSD of the crystal -> `results/overview/bioactive-cost/` (`bioactive_cost_ranked.png`,
`bioactive_cost_cumulative.png` (dG_dist only), `bioactive_cost.csv`). Crystals cluster_cut
cannot place (cut through a polar H; symmetric assignments in different clusters) are left out of the
figures but keep their RMSD-based values in the CSV.

## KL and Jensen-Shannon on ONE joint clustering

`kl_joint.py` (`bash kl_joint_submit.sh`: one 30-min task per compound, then the figures): barrier cuts
(defaults) on the water and complex frames together, the longer run subsampled to the shorter one's frame
count (equal weight: the clustered ensemble is m = (pi^w + pi^p)/2); every frame of both runs then labelled
through the joint cut tree. KL(p||w), KL(w||p) with a 0.5-frame pseudocount per cluster (clusters only one
ensemble visits; `*` in the figure where they hold > 1 % of one ensemble) and JS = (KL(w||m) + KL(p||m))/2,
each half in [0, ln 2] -> `results/overview/kl-divergence/kl_joint_cut.png`, `js_joint_halves_cut.png`,
`js_joint_ranked_cut.png` (CSVs in `data/`). Folate: one joint clustering per complex.
