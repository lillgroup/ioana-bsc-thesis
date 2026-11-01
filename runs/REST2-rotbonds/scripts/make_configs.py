"""Water and complex configs for every compound of compounds.csv, copied from the earlier
small-molecule v4 runs so that only the molecule changes:

  <name>-water/   build/build-top.config, build/scaler.config, REST2.config
                  = REST2-conformers-v4/4ibi-water (6 replicas, tau 0..0.4, 30000 exchanges = 60 ns)
  <name>-complex/ build/build-top.config, build/scaler.config, REST2.config
                  = REST2-complex-xtalwater/4ibi (crystal waters, ff14SB, PROPKA pH 7.4, TIP3P, 1.2 nm
                  cubic box; ligand-only REST2, 6 replicas, tau 0..0.4, 20000 exchanges = 40 ns)
The complex's ligand residue number is not known before the build: scaler.config and REST2.config
carry the mask ":LIGRES", which scripts/build.sbatch replaces with the built residue number.

  python -I scripts/make_configs.py <root> [name ...]    (names: only those rows; default all -- note that
  rewriting a built complex's configs resets its mask to :LIGRES)
"""
import csv
import sys
from pathlib import Path

WATER_BUILD = """# {name} ({ligand}) in water, v4: as REST2-conformers-v4/4ibi-water, with this compound's package.
solute:
  kind: ligand
  parameters: ../../packages/{name}
protonation:
  ligand_method: openbabel
  ph: 7.4
"""
WATER_SCALER = """method: REST2
schedule:
  kind: linear
  n_states: 6
  tau_min: 0.0
  tau_max: 0.4
"""
WATER_REST2 = """# {name} ({ligand}) in water: as REST2-conformers-v4/<n>-water, 30000 exchanges x 1000 steps x 2 fs = 60 ns per state.
protocol: REST2
solvent: explicit

rest2:
  number_of_replicas: 6
  tau_max: 0.4
  exchange_interval_steps: 1000
  number_of_exchanges: 30000
  equilibration_steps: 50000

stages:
  minimization_iterations: 5000
  restrained_nvt_steps: 50000
  restrained_npt_steps: 50000
  unrestrained_npt_steps: 100000

reporting:
  crd_printout_solute: 1000
  crd_printout_whole: 50000
  info_printout: 5000
  checkpoint_printout: 50000
"""
COMPLEX_BUILD = """# {pdb} {protein} chain {chain} + {ligand} (crystal pose, heavy atoms; H from the package) + crystal waters.
# Ligand parameters = the package of the water run (packages/{name}), so both runs share one Hamiltonian.
# Build settings as REST2-complex-xtalwater.
solute:
  kind: complex
ligands:
  - resname: {ligand}
    parameter: ../../packages/{name}
input:
  missing_atoms: add
protonation:
  method: propka
  ph: 7.4
forcefield:
  protein: ff14SB
solvent:
  model: TIP3P
  padding_nm: 1.2
  box_shape: cube
  cutoff_nm: 0.9
constraints:
  type: HBonds
  rigid_water: true
"""
COMPLEX_SCALER = """# Selective REST2: ligand only (residue :LIGRES = {ligand}, set by build.sbatch), 6 states, tau 0..0.4 linear (as the water run).
method: REST2
schedule: {{kind: linear, n_states: 6, tau_min: 0.0, tau_max: 0.4}}
ligand_scaling_dict:
  {ligand}:
    mask: ":LIGRES"
    torsion_exclusions: auto
"""
COMPLEX_REST2 = """# {name} complex ({ligand}): selective REST2, ligand only. 20000 exchanges x 1000 steps x 2 fs = 40 ns per state.
# As REST2-complex-xtalwater/4ibi; identical to the water run's REST2.config except the length and the
# ligand-only hot region.
protocol: REST2
solvent: explicit

rest2:
  number_of_replicas: 6
  tau_max: 0.4
  exchange_interval_steps: 1000
  number_of_exchanges: 20000
  equilibration_steps: 50000
  state_trajectory: true
  rem_log: true
  neighbour_acceptance_report: true
  ligand_scaling_dict:
    {ligand}: {{mask: ":LIGRES"}}

stages:
  minimization_iterations: 5000
  restrained_nvt_steps: 50000
  restrained_npt_steps: 50000
  unrestrained_npt_steps: 100000

reporting:
  crd_printout_solute: 1000
  crd_printout_whole: 50000
  info_printout: 5000
  checkpoint_printout: 50000
"""

if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve()
    names = set(sys.argv[2:])
    for r in csv.DictReader(open(root / "compounds.csv")):
        if names and r["name"] not in names:
            continue
        files = {
            f"{r['name']}-water/build/build-top.config": WATER_BUILD,
            f"{r['name']}-water/build/scaler.config": WATER_SCALER,
            f"{r['name']}-water/REST2.config": WATER_REST2,
            f"{r['name']}-complex/build/build-top.config": COMPLEX_BUILD,
            f"{r['name']}-complex/build/scaler.config": COMPLEX_SCALER,
            f"{r['name']}-complex/REST2.config": COMPLEX_REST2,
        }
        for path, text in files.items():
            (root / path).parent.mkdir(parents=True, exist_ok=True)
            (root / path).write_text(text.format(**r))
        print(f"{r['name']}: {len(files)} configs")
