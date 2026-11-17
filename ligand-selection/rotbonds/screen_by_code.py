"""The structure screen of screen_candidates.py for candidates whose PDBbind files are not at hand (on
sciCORE): the ligand is found in the deposited PDB entry by its residue code (refined_set_properties.csv's
ligname), first instance with altloc ' '/'A', instead of by PDBbind's coordinates. Every other check is the
same: chain breaks in the kept chains, hetero groups in that chain (additives / ions vs other), other
polymer chains within 5 A (short ones kept as bound peptides, long ones = packing), crystal waters.
  python -I screen_by_code.py <candidates.txt: pdb rb mw ligname n_res uniprot> <properties.csv> <crystal-dir> <out.csv>
"""
import collections
import csv
import math
import sys

CAND, PROPS, CRYSTAL, OUT = sys.argv[1:5]
ADDITIVES = {"GOL","EDO","SO4","PO4","DMS","NO3","PG4","PEG","ACY","ACT","FMT","MPD","BME","EPE","TRS","IMD","MES","PGE","1PE","P6G","CIT","TLA","SCN","IOD","BU3","DTT","MLI","NH4","EOH","IPA","PE4"}
IONS = {"NA","CL","K","MG","CA","ZN","BR","I"}
xyz = lambda l: (float(l[30:38]), float(l[38:46]), float(l[46:54]))
props = {r["pdb"]: r for r in csv.DictReader(open(PROPS))}
rows = []
for line in open(CAND):
    code, rb, mw, lig, n_res, uni = line.split()[:6]
    P = code.upper()
    try:
        recs = [l for l in open(f"{CRYSTAL}/{P}.pdb") if l.startswith(("ATOM", "HETATM")) and l[16] in " A"
                and l[76:78].strip() not in ("H", "D")]
    except FileNotFoundError:
        rows.append(dict(pdb=code, rb=rb, note="no PDB-format file", ok=False)); continue
    groups = collections.OrderedDict()
    for l in recs:
        if l.startswith("HETATM"):
            groups.setdefault((l[17:20].strip(), l[21], l[22:27]), []).append(l)
    hits = [k for k in groups if k[0] == lig]
    if not hits:
        rows.append(dict(pdb=code, rb=rb, lig=lig, note=f"ligand {lig} not found as a HETATM residue", ok=False)); continue
    best = hits[0]
    ch = best[1]
    lc = [xyz(l) for l in groups[best]]
    near = lambda ls, cut: min(math.dist(xyz(l), q) for l in ls for q in lc) < cut
    polymer = collections.defaultdict(list)
    for l in recs:
        if l.startswith("ATOM"):
            polymer[l[21]].append(l)
    contact = [k for k, ls in polymer.items() if k != ch and near(ls, 5.0)]
    nres = lambda k: len({l[22:27] for l in polymer[k]})
    other_chains = [k for k in contact if nres(k) < 30]
    packing = [k for k in contact if nres(k) >= 30]
    keep = [ch] + other_chains
    breaks = []
    for k in keep:
        res = collections.OrderedDict()
        for l in polymer[k]:
            res.setdefault(l[22:27], {})[l[12:16].strip()] = xyz(l)
        ids = list(res)
        for a, b in zip(ids, ids[1:]):
            if "C" in res[a] and "N" in res[b] and math.dist(res[a]["C"], res[b]["N"]) > 2.0:
                breaks.append(f"{k}{a.strip()}-{b.strip()}")
    hets, odd = [], []
    for (rn, c2, ri), ls in groups.items():
        if c2 not in keep or (rn, c2, ri) == best or rn == "HOH":
            continue
        dmin = min(math.dist(xyz(l), q) for l in ls for q in lc)
        tag = f"{rn}{ri.strip()}:{dmin:.0f}"
        if rn == lig:
            continue                       # other copies of the ligand in the kept chain: not kept anyway
        (hets if rn in ADDITIVES or rn in IONS else odd).append(tag)
    nw = sum(1 for l in recs if l[17:20] == "HOH" and l[21] in keep)
    p = props.get(code, {})
    rows.append(dict(pdb=code, rb=rb, mw=mw, lig=lig, chain="".join(keep), n_res=n_res, uniprot=uni,
                     protein=p.get("protein", "")[:28], res_A=p.get("resolution", ""), affinity=p.get("affinity", ""),
                     waters=nw, breaks=" ".join(breaks), additives_ions=" ".join(hets), other_het=" ".join(odd),
                     packing="".join(packing), ok=not breaks and not odd, note=""))
keys = ["pdb", "rb", "mw", "lig", "chain", "n_res", "uniprot", "protein", "res_A", "affinity", "waters", "breaks",
        "additives_ions", "other_het", "packing", "ok", "note"]
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in keys})
for r in rows:
    print(r["pdb"], r["rb"], r.get("lig", ""), r.get("n_res", ""), r.get("protein", "")[:24], r.get("res_A", ""),
          "OK" if r["ok"] else "--", r.get("breaks", "")[:30], r.get("other_het", "")[:40], r.get("note", ""))
