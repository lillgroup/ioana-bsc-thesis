#!/usr/bin/env python
"""
Render a trajectory movie (MP4) from a folder made by export_pymol.py, without opening PyMOL:
loads view.pml, renders one PNG per state headless, then ffmpeg joins them with a time stamp.

  python render_movie.py <export dir> [--out movie.mp4] [--fps 15] [--size 1280x720]
                         [--water] [--ray] [--every N]

--water keeps the waters (hidden by default), --ray ray-traces each frame (nicer, ~10x slower),
--every N renders only every N-th state (quick previews).
"""
import argparse
import os
import shutil
import subprocess
import tempfile

from pymol import cmd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("folder")
    ap.add_argument("--out", default=None)
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--size", default="1280x720")
    ap.add_argument("--water", action="store_true")
    ap.add_argument("--ray", action="store_true")
    ap.add_argument("--every", type=int, default=1)
    args = ap.parse_args()
    folder = os.path.abspath(args.folder)
    out = args.out or f"{folder}/{os.path.basename(folder)}.mp4"
    w, h = map(int, args.size.split("x"))

    os.chdir(folder)
    cmd.do("@view.pml")
    if not args.water:
        cmd.hide("spheres", "resn HOH")
    cmd.hide("labels", "contacts")
    cmd.set("ray_opaque_background", 1)
    cmd.set("max_threads", os.cpu_count() or 4)
    cmd.viewport(w, h)
    focus = "lig or cof" if cmd.count_atoms("cof") else "lig"
    cmd.orient(focus)
    cmd.zoom(focus, 2)
    cmd.clip("slab", 60)
    dt = None
    for line in open("view.pml"):
        if "state k =" in line:
            dt = float(line.split("state k =")[1].split("ns")[0])
            break

    n = cmd.count_states("md")
    tmp = tempfile.mkdtemp(prefix="frames_", dir=os.path.dirname(out))
    states = list(range(1, n + 1, args.every))
    for k, s in enumerate(states):
        cmd.frame(s)
        cmd.set("state", s)
        cmd.png(f"{tmp}/f{k:05d}.png", width=w, height=h, ray=1 if args.ray else 0, quiet=1)
        if k % 50 == 0:
            print(f"frame {k + 1}/{len(states)}", flush=True)

    step = (dt or 0.1) * args.every
    stamp = (f"drawtext=text='%{{eif\\:n*{step * 10:g}/10\\:d}}.%{{eif\\:mod(n*{step * 10:g}\\,10)\\:d}} ns':"
             f"x=20:y=h-50:fontsize=32:fontcolor=black")
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-framerate", str(args.fps), "-i", f"{tmp}/f%05d.png",
                    "-vf", f"{stamp},pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p",
                    "-crf", "18", "-movflags", "+faststart", out], check=True)
    shutil.rmtree(tmp)
    print(f"wrote {out}: {len(states)} frames, {len(states) / args.fps:.0f} s, {step:g} ns per frame")


if __name__ == "__main__":
    main()
