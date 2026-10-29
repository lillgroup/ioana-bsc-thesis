#!/usr/bin/env python
"""
Did the ligand stay in the pocket, and in its crystal pose? For every complex REST2 run, state 0,
per frame, four measures:

  RMSD    ligand heavy-atom RMSD to the crystal pose in the POCKET frame (A): each frame is
          superposed on the reference by its pocket C-alphas only (Kabsch), the same transform is
          applied to the ligand, and the ligand is NOT fitted again -- so the RMSD measures how far
          the ligand moved relative to the binding site. No symmetry correction (flipped
          carboxylate O's / ring flips add a little).
  COM     distance between the ligand centre of mass and the pocket centre of mass (A), computed
          within each frame (no alignment). Constant while bound even if the ligand turns; it
          detects leaving, not reorienting (a ligand on the pocket rim can look similar).
  contacts  all contacts: number of ligand heavy atom - receptor heavy atom pairs < 4 A
          (receptor = protein + cofactor, i.e. NADP+ in 3QL3). Zero = not touching the protein.
  Q       native contacts: of the pairs < 4 A in the crystal pose, the fraction still present
          (< 4.5 A, 0.5 A slack for thermal motion). ~1 crystal binding mode, ~0.5 partly
          rearranged, ~0 crystal interactions lost (even if still on the protein).

Reference = build/built.solute.pdb: the crystal complex as prepared for the simulation (same
coordinates as build/built.pdb without the water), so protein and ligand are already where the
crystal puts them. Pocket = protein residues with a heavy atom within 5 A of the ligand there.
Distances use the minimum image of the cubic box. Thin line = every analysed frame, thick = 1 ns
running mean.

Writes binding_folate.png, binding_small_molecules.png (one row per run) and binding_summary.csv
(whole run and last 10 %, plus the closest contact for a dissociation check).
Needs numpy, scipy, matplotlib (bsc-thesis/venv).
  python check_bound.py [--thesis DIR] [--out DIR] [--stride N]
"""
import argparse
import csv
import glob
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

HOME = os.path.expanduser("~")
RUNS = ["FOL-REST2-v4/4LRH", "FOL-REST2-v4/3QL3", "FOL-REST2-mace/4LRH", "FOL-REST2-mace/3QL3",
        "FOL-complex-REST2/4LRH", "FOL-complex-REST2/3QL3",
        "REST2-conformers-v4/4ibi-complex", "REST2-conformers-v4/4ibk-complex",
        "REST2-conformers-v4/4j7d-complex",
        "REST2-conformers-complex/4ibi", "REST2-conformers-complex/4ibk",
        "REST2-conformers-complex/4j7d", "REST2-conformers-complex/5mwh",
        "REST2-conformers-complex/5ov8"]
PROTEIN = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS HID HIE HIP ILE LEU LYS MET PHE PRO SER THR "
              "TRP TYR VAL CYX ACE NME".split())
POCKET_CUT, CONTACT_CUT, NATIVE_SLACK, UNBOUND_CUT = 5.0, 4.0, 0.5, 5.0
MASS = {"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974, "F": 18.998, "CL": 35.45,
        "BR": 79.904, "I": 126.90, "H": 1.008}
LABELS = {"FOL-REST2-v4": "v4 (Sage)", "FOL-REST2-mace": "MACE bespoke", "FOL-complex-REST2": "original",
          "REST2-conformers-v4": "v4", "REST2-conformers-complex": "original"}


def read_pdb(path):
    atoms, xyz = [], []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            el = line[76:78].strip() or line[12:16].strip()[0]
            atoms.append({"name": line[12:16].strip(), "res": line[17:20].strip(), "el": el.upper(),
                          "resid": (line[21], line[22:27]), "heavy": el.upper() != "H"})
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return atoms, np.array(xyz)


def min_image(d, box):
    return d - box * np.round(d / box)


def kabsch(P, Q):
    """Rotation R minimising |P R - Q| for centred P, Q."""
    U, _, Vt = np.linalg.svd(P.T @ Q)
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    return U @ D @ Vt


def analyse(run_dir, stride):
    ligand = os.path.splitext(os.path.basename(glob.glob(f"{run_dir}/build/*.sdf")[0]))[0]
    atoms, ref = read_pdb(f"{run_dir}/build/built.solute.pdb")
    lig = np.array([i for i, a in enumerate(atoms) if a["res"] == ligand and a["heavy"]])
    rec = np.array([i for i, a in enumerate(atoms) if a["res"] != ligand and a["heavy"]])
    box_ref = float(open(f"{run_dir}/build/built.solute.pdb").read().split("CRYST1")[1].split()[0])

    d_ref = np.linalg.norm(min_image(ref[lig][:, None] - ref[rec][None], box_ref), axis=-1)
    li, ri = np.nonzero(d_ref < CONTACT_CUT)  # native pairs
    near = {atoms[rec[j]]["resid"] for j in np.nonzero((d_ref < POCKET_CUT).any(0))[0]
            if atoms[rec[j]]["res"] in PROTEIN}
    pocket_ca = np.array([i for i, a in enumerate(atoms) if a["name"] == "CA" and a["resid"] in near])
    pocket = np.array([i for i, a in enumerate(atoms) if a["heavy"] and a["resid"] in near
                       and a["res"] in PROTEIN])
    w_lig = np.array([MASS.get(atoms[i]["el"], 12.0) for i in lig])
    w_poc = np.array([MASS.get(atoms[i]["el"], 12.0) for i in pocket])
    com = lambda x, w: (x * w[:, None]).sum(0) / w.sum()

    f = netcdf_file(f"{run_dir}/REST2-run1/solute_state0_prod1.nc", "r", mmap=True)
    xyz = np.array(f.variables["coordinates"][::stride], dtype=float)
    boxes = np.array(f.variables["cell_lengths"][::stride], dtype=float)[:, 0]
    t_ns = np.array(f.variables["time"][::stride], dtype=float) / 1000.0
    f.close()

    ref_ca_c = ref[pocket_ca].mean(0)
    ref_lig = ref[lig] - ref_ca_c
    out = {k: np.empty(len(xyz)) for k in ("rmsd", "com", "contacts", "Q", "d_min")}
    for n, (x, box) in enumerate(zip(xyz, boxes)):
        # pocket made whole around its first C-alpha, ligand brought to the nearest image of the pocket
        anchor = x[pocket_ca[0]]
        ca = anchor + min_image(x[pocket_ca] - anchor, box)
        poc = anchor + min_image(x[pocket] - anchor, box)
        ca_c = ca.mean(0)
        lg = x[lig]
        lg = lg[0] + min_image(lg - lg[0], box)
        lg = lg + (ca_c + min_image(lg.mean(0) - ca_c, box) - lg.mean(0))
        d = np.linalg.norm(min_image(lg[:, None] - x[rec][None], box), axis=-1)
        out["d_min"][n] = d.min()
        out["contacts"][n] = (d < CONTACT_CUT).sum()
        out["Q"][n] = np.mean(d[li, ri] < CONTACT_CUT + NATIVE_SLACK) if len(li) else np.nan
        out["com"][n] = np.linalg.norm(com(lg, w_lig) - com(poc, w_poc))
        R = kabsch(ca - ca_c, ref[pocket_ca] - ref_ca_c)  # pocket fit only, ligand not refitted
        out["rmsd"][n] = np.sqrt((((lg - ca_c) @ R - ref_lig) ** 2).sum(1).mean())
    com0 = np.linalg.norm(com(ref[lig], w_lig) - com(ref[pocket], w_poc))
    native = {"pairs": len(li), "contacts0": int((d_ref < CONTACT_CUT).sum()), "com0": com0}
    return ligand, t_ns, out, len(near), native


def running_mean(y, n):
    if n < 2 or len(y) < n:
        return y
    k = np.ones(n) / n
    return np.convolve(np.pad(y, (n // 2, n - 1 - n // 2), mode="edge"), k, mode="valid")


def run_label(run, ligand):
    group, system = run.split("/")
    system = system.replace("-complex", "")
    name = "folate" if ligand == "FOL" else ligand
    return f"{name} · {system}\n{LABELS.get(group, group)}"


COLS = (("rmsd", "ligand RMSD, pocket frame (Å)"), ("com", "ligand–pocket COM distance (Å)"),
        ("contacts", "all contacts (< 4 Å)"), ("Q", "native contacts Q"))


def plot_group(results, title, path):
    fig, axes = plt.subplots(len(results), len(COLS), figsize=(17, 1.75 * len(results) + 1.2),
                             sharex=True, sharey="col", squeeze=False)
    for r, (row, t, m, native) in enumerate(results):
        win = max(1, int(round(1.0 / (t[1] - t[0]))))  # 1 ns
        for c, (key, label) in enumerate(COLS):
            ax = axes[r, c]
            ax.plot(t, m[key], lw=0.4, color="#86b6ef")
            ax.plot(t, running_mean(m[key], win), lw=1.6, color="#1c5cab")
            if key == "com":
                ax.axhline(native["com0"], color="#555", lw=0.8, ls=":")
            if key == "contacts":
                ax.axhline(native["contacts0"], color="#555", lw=0.8, ls=":")
            if key == "Q":
                ax.set_ylim(0, 1.02)
            ax.tick_params(labelsize=8)
            ax.grid(axis="y", color="#e6e6e6", lw=0.6)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            if r == 0:
                ax.set_title(label, fontsize=11)
            if r == len(results) - 1:
                ax.set_xlabel("time (ns)", fontsize=9)
        axes[r, 0].set_ylabel(f"{run_label(row['run'], row['ligand'])}\n{row['verdict']}", fontsize=8.5,
                              rotation=0, ha="right", va="center", labelpad=12)
    for ax in axes[:, 0]:
        ax.set_ylim(bottom=0)
    for ax in axes[:, 1]:
        ax.set_ylim(bottom=0)
    for ax in axes[:, 2]:
        ax.set_ylim(bottom=0)
    fig.suptitle(title, fontsize=13, y=0.995)
    fig.text(0.5, 0.003, "state 0 · thin: every analysed frame, thick: 1 ns running mean · dotted: crystal value · "
             "same y-scale down each column", ha="center", fontsize=9, color="#555")
    fig.tight_layout(rect=(0, 0.015, 1, 0.98))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--thesis", default=f"{HOME}/thesis-sbrinz-data/thesis")
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--stride", type=int, default=10, help="analyse every N-th saved frame")
    args = ap.parse_args()

    results = []
    for run in RUNS:
        d = f"{args.thesis}/{run}"
        if not os.path.exists(f"{d}/REST2-run1/solute_state0_prod1.nc"):
            print(f"skip {run}: no trajectory")
            continue
        ligand, t, m, npocket, native = analyse(d, args.stride)
        tail = t >= 0.9 * t[-1]
        unbound = m["d_min"] > UNBOUND_CUT
        row = {"run": run, "ligand": ligand, "ns": round(t[-1], 1), "pocket_residues": npocket,
               "native_pairs": native["pairs"], "com_crystal": round(native["com0"], 2),
               "frac_unbound": round(unbound.mean(), 3), "max_closest_contact": round(m["d_min"].max(), 2)}
        for k in ("rmsd", "com", "contacts", "Q"):
            row[f"{k}_mean"] = round(float(np.nanmean(m[k])), 2)
            row[f"{k}_last10pct"] = round(float(np.nanmean(m[k][tail])), 2)
        row["verdict"] = ("DISSOCIATED" if unbound.any() else
                          "bound, crystal pose lost" if row["Q_last10pct"] < 0.3 else
                          "bound, rearranged" if row["Q_last10pct"] < 0.6 or row["rmsd_last10pct"] > 3 else
                          "bound, crystal pose")
        results.append((row, t, m, native))
        print(f"{run:36s} {ligand}  RMSD {row['rmsd_last10pct']:5.2f} A  COM {row['com_last10pct']:5.2f} A "
              f"(crystal {row['com_crystal']:5.2f})  contacts {row['contacts_last10pct']:5.1f}  "
              f"Q {row['Q_last10pct']:.2f}  -> {row['verdict']}")

    with open(f"{args.out}/binding_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0][0]))
        w.writeheader()
        w.writerows(r[0] for r in results)

    fol = [r for r in results if r[0]["ligand"] == "FOL"]
    small = [r for r in results if r[0]["ligand"] != "FOL"]
    plot_group(fol, "Folate in folate receptor α (4LRH) and DHFR (3QL3)", f"{args.out}/binding_folate.png")
    plot_group(small, "Small molecules in their proteins", f"{args.out}/binding_small_molecules.png")
    print(f"wrote {args.out}/binding_folate.png, binding_small_molecules.png and binding_summary.csv")



if __name__ == "__main__":
    main()
