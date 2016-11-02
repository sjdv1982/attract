#fills a PDB B-factor (temperature) column from per-residue numerical data
#input: PDB file, numerical data
#column 1-6 of the numerical data is the full resid (chainID + 5 resid columns; chainID _ is interpreted as " ")
#the rest of the line is interpreted as a numerical value (float), V
# B-factor = V / max(V) * 100
# i.e. the largest numerical value (max(V)) gets a B-factor of 100
# missing numerical data is interpreted as zero

import sys, rmsdlib
pdb = rmsdlib.read_pdb(sys.argv[1])
resids = [res[0].resid for res in pdb.residues()]
values = {}
for l in open(sys.argv[2]):
    resid = l[:6].replace("_"," ")
    assert resid in resids, resid
    v = float(l[6:])
    values[resid] = v
maxv = max(values.values())

for res in pdb.residues():
    resid = res[0].resid
    v = values.get(resid, 0.0)
    b = v/maxv * 100.0
    for atom in res:
        l = atom.line
        ll = l[:60] + "%6.2f" % b + l[66:]
        sys.stdout.write(ll)