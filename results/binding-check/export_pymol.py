#!/usr/bin/env python
"""
Export one complex REST2 run (state 0) as PyMOL-ready files: start.pdb (the prepared crystal
complex), md.pdb (multi-state movie) and view.pml (representations, pocket, contacts, pose clouds).

Frames come from REST2-run1/whole_state0_prod1.nc (full system incl. water, every 100 ps).
Every frame is made whole (walking the atom list with the minimum image of the cubic box) and
superposed on the pocket C-alphas, so in the movie the protein stays still and only the ligand's
own motion is left. Protein H are dropped. Water: only the --waters closest to the ligand in each
frame are kept (as fixed "slots" W1..Wn, ordered by distance, so a slot is not the same molecule
from frame to frame -- fine for seeing the hydration shell, not for following one water).

  python export_pymol.py [--run DIR] [--ligand FOL] [--waters 200] [--out DIR]
then in PyMOL:  cd <out>   and   @view.pml
"""
import argparse
import os
import sys

import numpy as np
from scipy.io import netcdf_file

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_bound import read_pdb, min_image, kabsch, PROTEIN, POCKET_CUT

HOME = os.path.expanduser("~")
KEY_RESIDUES = {"3QL3": "(resn ARG and resi 57) or (resn ASP and resi 27)"}


def pdb_lines(path):
    return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def fmt(line, xyz, serial):
    return f"{line[:6]}{serial % 100000:5d}{line[11:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", default=f"{HOME}/thesis-sbrinz-data/thesis/FOL-REST2-v4/3QL3")
    ap.add_argument("--ligand", default="FOL")
    ap.add_argument("--waters", type=int, default=200)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    run = args.run.rstrip("/")
    name = "_".join(run.split("/")[-2:])
    out = args.out or f"{HOME}/results_bsc_thesis/pymol/{name}"
    os.makedirs(out, exist_ok=True)

    top = f"{run}/build/built.pdb"
    atoms, ref = read_pdb(top)
    lines = pdb_lines(top)
    box0 = float(open(top).read().split("CRYST1")[1].split()[0])
    lig = np.array([i for i, a in enumerate(atoms) if a["res"] == args.ligand])
    lig_heavy = np.array([i for i in lig if atoms[i]["heavy"]])
    water_o = np.array([i for i, a in enumerate(atoms) if a["res"] == "HOH" and a["name"] == "O"])
    keep = np.array([i for i, a in enumerate(atoms) if a["res"] not in ("HOH", "NA", "CL", "K")
                     and (a["heavy"] or a["res"] == args.ligand)])
    rec_heavy = np.array([i for i in keep if atoms[i]["heavy"] and atoms[i]["res"] != args.ligand])
    d0 = np.linalg.norm(min_image(ref[lig_heavy][:, None] - ref[rec_heavy][None], box0), axis=-1)
    near = {atoms[rec_heavy[j]]["resid"] for j in np.nonzero((d0 < POCKET_CUT).any(0))[0]
            if atoms[rec_heavy[j]]["res"] in PROTEIN}
    ca = np.array([i for i, a in enumerate(atoms) if a["name"] == "CA" and a["resid"] in near])
    ref_c = ref[ca].mean(0)
    kpos = {a: k for k, a in enumerate(keep)}
    ca_k = np.array([kpos[i] for i in ca])
    het_groups = [np.array([kpos[i] for i in keep if atoms[i]["resid"] == r and atoms[i]["res"] not in PROTEIN])
                  for r in dict.fromkeys(atoms[i]["resid"] for i in keep if atoms[i]["res"] not in PROTEIN)]

    with open(f"{out}/start.pdb", "w") as fh:
        for n, i in enumerate(keep, 1):
            fh.write(fmt(lines[i], ref[i], n))
        fh.write("END\n")

    f = netcdf_file(f"{run}/REST2-run1/whole_state0_prod1.nc", "r", mmap=False)
    xyz_all = f.variables["coordinates"]
    boxes = f.variables["cell_lengths"][:, 0].astype(float)
    t_ns = f.variables["time"][:].astype(float) / 1000.0
    water_tmpl = [lines[water_o[0]], lines[water_o[0] + 1], lines[water_o[0] + 2]]
    with open(f"{out}/md.pdb", "w") as fh:
        for s in range(len(t_ns)):
            x = np.array(xyz_all[s], dtype=float)
            box = boxes[s]
            # make whole by walking the atom list (consecutive atoms within a molecule are bonded or
            # neighbours, < half a box apart), then put each ligand/cofactor at the image nearest
            # the pocket -- the jump from the protein's last atom to them can exceed half a box
            xk = x[keep].copy()
            xk[1:] = min_image(np.diff(xk, axis=0), box)
            xk = np.cumsum(xk, axis=0)
            pocket_c = xk[ca_k].mean(0)
            for g in het_groups:
                com = xk[g].mean(0)
                xk[g] += pocket_c + min_image(com - pocket_c, box) - com
            pos = dict(zip(keep, xk))
            lig_now = np.array([pos[i] for i in lig_heavy])
            # waters: nearest image of each water O to the ligand, keep the closest N (with their H)
            wo = x[water_o]
            lc = lig_now.mean(0)
            wo = lc + min_image(wo - lc, box)
            dist = np.linalg.norm(wo[:, None] - lig_now[None], axis=-1).min(1)
            pick = np.argsort(dist)[: args.waters]
            wat = []
            for w in pick:
                o = water_o[w]
                shift = wo[w] - x[o]
                wat.append([x[o] + shift, x[o + 1] + shift, x[o + 2] + shift])
            # superpose on the pocket C-alphas
            P = np.array([pos[i] for i in ca])
            pc = P.mean(0)
            R = kabsch(P - pc, ref[ca] - ref_c)
            move = lambda v: (v - pc) @ R + ref_c
            fh.write(f"MODEL     {s + 1:4d}\nREMARK time {t_ns[s]:.1f} ns\n")
            serial = 0
            for i in keep:
                serial += 1
                fh.write(fmt(lines[i], move(pos[i]), serial))
            for k, w in enumerate(wat):
                for line, v in zip(water_tmpl, w):
                    serial += 1
                    l = line[:17] + "HOH W" + f"{k + 1:4d}" + line[26:]
                    fh.write(fmt(l, move(v), serial))
            fh.write("ENDMDL\n")
    f.close()
    dt = t_ns[1] - t_ns[0]

    key = next((v for k, v in KEY_RESIDUES.items() if k in run), None)
    if key is None:
        # no curated list: residues with an N/O within 3.5 A of a ligand N/O in the starting pose
        polar = lambda i: atoms[i]["name"][0] in "NO"
        lp = [i for i in lig_heavy if polar(i)]
        hits = {}
        for j in rec_heavy:
            if polar(j) and atoms[j]["res"] in PROTEIN and \
                    np.linalg.norm(min_image(ref[lp] - ref[j], box0), axis=1).min() < 3.5:
                hits[atoms[j]["resid"]] = atoms[j]["res"]
        key = " or ".join(f"(resn {res} and resi {rid[1].strip()})" for rid, res in hits.items()) or None
        print("key residues (polar contacts in the starting pose):",
              ", ".join(f"{res}{rid[1].strip()}" for rid, res in hits.items()))
    late0 = int(0.8 * len(t_ns)) + 1
    pml = f"""# {name}: {args.ligand}, REST2 state 0, {t_ns[-1]:.0f} ns, {len(t_ns)} states (state k = {dt:.1f} ns x k)
# start  = prepared crystal complex (reference pose, white sticks)
# md     = the movie, pocket-aligned, {args.waters} closest waters per frame (red spheres)
# cloud_early / cloud_late = ligand poses overlaid, first vs last 20 % of the run
load start.pdb, start
load md.pdb, md
bg_color white
hide everything
set cartoon_transparency, 0.55
show cartoon, md and polymer.protein
color grey80, md and polymer.protein and elem C
select lig, md and resn {args.ligand}
select cof, md and organic and not resn {args.ligand}
select pocket, byres (md and polymer.protein) within 4.5 of (start and resn {args.ligand})
show sticks, (pocket and not name N+C+O) or lig or cof
util.cbag lig
util.cbay cof
color grey60, pocket and elem C
hide everything, start
show sticks, start and resn {args.ligand} and not hydro
color white, start and resn {args.ligand} and elem C
set stick_transparency, 0.5, start
show spheres, md and resn HOH and elem O
set sphere_scale, 0.22, md and resn HOH
color red, md and resn HOH and elem O
hide everything, lig and hydro and (elem H and neighbor elem C)
{f'''select key, md and ({key})
show sticks, key and not name N+C+O
color magenta, key and elem C
label key and name CA, "%s%s" % (resn, resi)''' if key else ''}
distance contacts, lig, (md and not lig and not hydro), mode=2, state=0
color yellow, contacts
hide labels, contacts
python
n = cmd.count_states("md")
step = max(1, n // 100)
early = range(1, int(0.2 * n) + 1, step)
late = range({late0}, n + 1, step)
for obj, states in (("cloud_early", early), ("cloud_late", late)):
    for k, s in enumerate(states, 1):
        cmd.create(obj, "md and resn {args.ligand} and not hydro", s, k)
    cmd.set("all_states", 1, obj)
    cmd.hide("everything", obj)
    cmd.show("lines", obj)
cmd.color("skyblue", "cloud_early")
cmd.color("orange", "cloud_late")
cmd.disable("cloud_early")
cmd.disable("cloud_late")
python end
orient lig
set movie_fps, 15
# play: the > button bottom right (or mplay). Clouds: click cloud_early / cloud_late in the panel
"""
    open(f"{out}/view.pml", "w").write(pml)
    print(f"wrote {out}/start.pdb, md.pdb ({len(t_ns)} states, {os.path.getsize(f'{out}/md.pdb') / 1e6:.0f} MB), view.pml")


if __name__ == "__main__":
    main()
