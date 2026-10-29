#!/usr/bin/env python
"""
Label every state-0 frame of each complex REST2 run as BOUND or UNBOUND from the ligand-pocket
centre-of-mass distance, so later analyses (torsions, clustering) can use bound frames only.

Why per frame: the state-0 trajectory is stitched from all replicas (an exchange every 2 ps), and a
replica whose ligand has left the pocket keeps passing through state 0. Its frames are mixed in
between bound frames, so whole time windows can't be cut; each frame is labelled on its own.

  COM distance = |COM(ligand heavy atoms) - COM(pocket heavy atoms)| within the frame (mass-
  weighted, minimum image). Pocket = protein residues within 5 A of the ligand in
  build/built.solute.pdb. Bound = COM distance <= --cutoff (default 6 A: the gap between the
  bound and the unbound peak, see com_histograms.png).

Writes, per run, <out>/<group>_<system>_com.csv (frame, time_ps, walker, com_A, bound), plus
com_histograms.png (one panel per run, cutoff marked) and com_fraction_bound.csv.
The walker column says which replica was in state 0 (from REST2-run1/exchange.csv).
Needs numpy, scipy, matplotlib (bsc-thesis/venv).
  python com_filter.py [--thesis DIR] [--out DIR] [--cutoff 6.0]
"""
import argparse
import csv
import glob
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.io import netcdf_file

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_bound import read_pdb, min_image, PROTEIN, POCKET_CUT, MASS, RUNS, run_label

HOME = os.path.expanduser("~")


def com_series(run_dir, chunk=2000):
    ligand = os.path.splitext(os.path.basename(glob.glob(f"{run_dir}/build/*.sdf")[0]))[0]
    atoms, ref = read_pdb(f"{run_dir}/build/built.solute.pdb")
    lig = np.array([i for i, a in enumerate(atoms) if a["res"] == ligand and a["heavy"]])
    rec = np.array([i for i, a in enumerate(atoms) if a["res"] != ligand and a["heavy"]])
    box0 = float(open(f"{run_dir}/build/built.solute.pdb").read().split("CRYST1")[1].split()[0])
    d = np.linalg.norm(min_image(ref[lig][:, None] - ref[rec][None], box0), axis=-1)
    near = {atoms[rec[j]]["resid"] for j in np.nonzero((d < POCKET_CUT).any(0))[0]
            if atoms[rec[j]]["res"] in PROTEIN}
    poc = np.array([i for i, a in enumerate(atoms) if a["heavy"] and a["resid"] in near and a["res"] in PROTEIN])
    wl = np.array([MASS.get(atoms[i]["el"], 12.0) for i in lig])
    wp = np.array([MASS.get(atoms[i]["el"], 12.0) for i in poc])

    f = netcdf_file(f"{run_dir}/REST2-run1/solute_state0_prod1.nc", "r", mmap=False)
    X, B = f.variables["coordinates"], f.variables["cell_lengths"]
    t = np.array(f.variables["time"][:], dtype=float)
    com = np.empty(len(t))
    for s in range(0, len(t), chunk):
        x = np.array(X[s:s + chunk], dtype=float)
        b = np.array(B[s:s + chunk], dtype=float)[:, 0]
        p = x[:, poc]
        p = p[:, :1] + min_image(p - p[:, :1], b[:, None, None])
        l = x[:, lig]
        l = l[:, :1] + min_image(l - l[:, :1], b[:, None, None])
        pc = (p * wp[None, :, None]).sum(1) / wp.sum()
        lc = (l * wl[None, :, None]).sum(1) / wl.sum()
        com[s:s + chunk] = np.linalg.norm(min_image(lc - pc, b[:, None]), axis=1)
    f.close()
    com0 = np.linalg.norm((ref[lig] * wl[:, None]).sum(0) / wl.sum() - (ref[poc] * wp[:, None]).sum(0) / wp.sum())

    walker = np.full(len(t), -1)
    exch = f"{run_dir}/REST2-run1/exchange.csv"
    if os.path.exists(exch):
        by_time = {}
        for r in csv.DictReader(open(exch)):
            if r["state"] == "0":
                by_time[round(float(r["time_ps"]), 1)] = int(r["walker"])
        walker = np.array([by_time.get(round(x, 1), -1) for x in t])
    return ligand, t, com, walker, com0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--thesis", default=f"{HOME}/thesis-sbrinz-data/thesis")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "com_filter"))
    ap.add_argument("--cutoff", type=float, default=6.0, help="bound if COM distance <= cutoff (A)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    results, summary = [], []
    for run in RUNS:
        d = f"{args.thesis}/{run}"
        if not os.path.exists(f"{d}/REST2-run1/solute_state0_prod1.nc"):
            continue
        ligand, t, com, walker, com0 = com_series(d)
        bound = com <= args.cutoff
        name = run.replace("/", "_")
        with open(f"{args.out}/{name}_com.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["frame", "time_ps", "walker", "com_A", "bound"])
            for i in range(len(t)):
                w.writerow([i, f"{t[i]:.1f}", walker[i], f"{com[i]:.3f}", int(bound[i])])
        unbound_walkers = sorted(set(walker[~bound]) - {-1})
        gap = np.mean((com > args.cutoff - 1) & (com < args.cutoff + 1))
        row = {"run": run, "ligand": ligand, "frames": len(t), "com_crystal": round(com0, 2),
               "frac_bound": round(bound.mean(), 4), "frac_within_1A_of_cutoff": round(gap, 4),
               "max_com": round(com.max(), 2), "unbound_frames_from_walkers": " ".join(map(str, unbound_walkers))}
        summary.append(row)
        results.append((run, ligand, com, com0))
        print(f"{run:36s} bound {bound.mean() * 100:6.2f} %   COM max {com.max():5.2f} A   "
              f"frames within +-1 A of cutoff {gap * 100:4.2f} %   unbound from walkers {unbound_walkers}")

    with open(f"{args.out}/com_fraction_bound.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)

    ncol = 4
    nrow = int(np.ceil(len(results) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(16, 2.6 * nrow + 0.8), sharex=True, squeeze=False)
    edges = np.arange(0, max(16.0, max(r[2].max() for r in results) + 0.5), 0.25)
    for ax, (run, ligand, com, com0) in zip(axes.flat, results):
        ax.hist(com[com <= args.cutoff], bins=edges, color="#2a78d6", label="bound")
        ax.hist(com[com > args.cutoff], bins=edges, color="#eb6834", label="unbound")
        ax.axvline(args.cutoff, color="#333", lw=1, ls="--")
        ax.axvline(com0, color="#333", lw=1, ls=":")
        ax.set_yscale("log")
        ax.set_title(f"{run_label(run, ligand).replace(chr(10), ' · ')}\n"
                     f"bound {np.mean(com <= args.cutoff) * 100:.1f} %", fontsize=9.5)
        ax.tick_params(labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    for ax in axes.flat[len(results):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("ligand–pocket COM distance (Å)", fontsize=9)
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(f"State-0 frames by ligand–pocket COM distance (log counts) · dashed: cutoff {args.cutoff:g} Å · "
                 "dotted: crystal", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(f"{args.out}/com_histograms.png", dpi=140)
    print(f"wrote {args.out}/ (per-run CSVs, com_fraction_bound.csv, com_histograms.png)")


if __name__ == "__main__":
    main()
