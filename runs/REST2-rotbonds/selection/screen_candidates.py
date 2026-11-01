"""Step 2 of the REST2-rotbonds selection: filter and screen the candidates for MD preparation.

  python -I selection/screen_candidates.py <refined-set> selection/refined_set_properties.csv <crystal-dir>

Filter (filter.py's criteria, with rotatable bonds counted on heavy atoms): MW 450-550, formal charge
0 in the PDBbind SDF (NB: the packages are re-protonated at pH 7.4 anyway), <= 40 heavy atoms, no
metal or cofactor in the pocket file, one protein chain, no B/Si/Se/P; 4-7 heavy-atom rotatable bonds.
Then per candidate: ligand chain (matching PDBbind's coordinates), internal chain breaks in the chains
kept, hetero groups in that chain (with distance to the ligand), other polymer chains within 5 A of
the ligand, crystal waters. Downloads the PDB entry from RCSB if missing."""
import math, os, sys, urllib.request, collections
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
R, PROPS, CRYSTAL = sys.argv[1], sys.argv[2], sys.argv[3]
ADDITIVES = {"GOL","EDO","SO4","PO4","DMS","NO3","PG4","PEG","ACY","ACT","FMT","MPD","BME","EPE","TRS","IMD","MES","PGE","1PE","P6G","CIT","TLA","SCN","IOD","BU3","DTT","MLI","NH4","EOH","IPA","PE4"}
IONS = {"NA","CL","K","MG","CA","ZN","BR","I"}
xyz = lambda l: (float(l[30:38]), float(l[38:46]), float(l[46:54]))
d = pd.read_csv(PROPS)
d = d[(d.mw.between(450, 550)) & (d.charge == 0) & (d.ha <= 40) & d.pocket_het.isna() & (d.n_chains == 1)
      & ~d.elements.str.contains("H|B|Si|Se|P")]
d = d[d.rb_heavy.isin([4,5,6,7])].sort_values(["rb_heavy","n_res"])
rows = []
for _, c in d.iterrows():
    code = c.pdb; P = code.upper(); f = f"{CRYSTAL}/{P}.pdb"
    if not os.path.exists(f):
        try: urllib.request.urlretrieve(f"https://files.rcsb.org/download/{P}.pdb", f)
        except Exception: rows.append(dict(pdb=code, note="no PDB-format file")); continue
    recs = [l for l in open(f) if l.startswith(("ATOM","HETATM")) and l[16] in " A" and l[76:78].strip() not in ("H","D")]
    m = Chem.MolFromMolFile(f"{R}/{code}/{code}_ligand.sdf", removeHs=False)
    pb = [tuple(p) for p, a in zip(m.GetConformer().GetPositions(), m.GetAtoms()) if a.GetAtomicNum() > 1]
    # ligand residue = the HETATM residue whose atoms coincide with PDBbind's
    groups = collections.defaultdict(list)
    for l in recs:
        if l.startswith("HETATM"): groups[(l[17:20].strip(), l[21], l[22:27])].append(l)
    def match(ls): return sum(min(math.dist(p, xyz(l)) for l in ls) < 0.1 for p in pb)
    best = max(groups, key=lambda k: match(groups[k])) if groups else None
    if best is None or match(groups[best]) < 0.9 * len(pb):
        rows.append(dict(pdb=code, note="ligand not found as one HETATM residue")); continue
    lig, ch = best[0], best[1]
    lc = [xyz(l) for l in groups[best]]
    near = lambda ls, cut: min(math.dist(xyz(l), q) for l in ls for q in lc) < cut
    polymer = collections.defaultdict(list)
    for l in recs:
        if l.startswith("ATOM"): polymer[l[21]].append(l)
    contact = [k for k, ls in polymer.items() if k != ch and near(ls, 5.0)]
    nres = lambda k: len({l[22:27] for l in polymer[k]})
    other_chains = [k for k in contact if nres(k) < 30]           # bound peptides (co-activator, ...)
    packing = [k for k in contact if nres(k) >= 30]
    keep = [ch] + other_chains
    breaks = []
    for k in keep:
        res = collections.OrderedDict()
        for l in polymer[k]:
            res.setdefault((l[22:26], l[26]), {})[l[12:16].strip()] = xyz(l)
        ks = list(res)
        breaks += [f"{k}{a[0].strip()}-{b[0].strip()}" for a, b in zip(ks, ks[1:])
                   if "C" in res[a] and "N" in res[b] and math.dist(res[a]["C"], res[b]["N"]) > 2.0]
    hets, odd = [], []
    for (rn, c2, ri), ls in groups.items():
        if c2 not in keep or (rn, c2, ri) == best or rn == "HOH": continue
        dmin = min(math.dist(xyz(l), q) for l in ls for q in lc)
        tag = f"{rn}{ri.strip()}:{dmin:.0f}"
        (hets if rn in ADDITIVES or rn in IONS else odd).append(tag)
    nw = sum(1 for l in recs if l[17:20] == "HOH" and l[21] in keep)
    rows.append(dict(pdb=code, rb=c.rb_heavy, lig=lig, chain="".join(keep), n_res=c.n_res,
                     protein=c.protein[:28], res_A=c.resolution, waters=nw, breaks=" ".join(breaks),
                     additives_ions=" ".join(hets), other_het=" ".join(odd),
                     packing="".join(packing), ok=not breaks and not odd))
out = pd.DataFrame(rows); out.to_csv("candidates_screened.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 40); pd.set_option("display.max_rows", 200)
out["protein"] = out.protein.str[:20]
for col in ("additives_ions", "other_het", "breaks"): out[col] = out[col].fillna("").str[:34]
print(out[["pdb","rb","lig","chain","packing","n_res","protein","res_A","waters","breaks","additives_ions","other_het","ok"]].to_string(index=False))
