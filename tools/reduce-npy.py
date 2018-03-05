import sys, os
import numpy as np
npy = np.load(sys.argv[1])
assert npy.ndim == 3, npy.ndim
pdb = sys.argv[2]
outpf = sys.argv[3]

dirname = os.path.dirname(os.path.abspath(__file__))

def read_forcefield(forcefieldfile):
    ff = {}
    aa = None
    for l in open(forcefieldfile):
        pound = l.find("#")
        if pound > -1: l = l[:pound]
        l = l.strip()
        if not len(l): continue
        ll = l.split()
        if len(ll) == 1:
            aa = ll[0]
            assert len(aa) <= 3, l
            ff[aa] = []
        else:
            assert aa is not None
            try:
                atomtype = int(ll[0])
            except ValueError:
                raise ValueError(l)
            atoms = ll[2:]
            charge = 0.0
            try:
                charge = float(atoms[-1])
                atoms = atoms[:-1]
            except ValueError:
                pass
            ff[aa].append( (int(ll[0]), ll[1], set(atoms), charge) )
    return ff
ff = read_forcefield(dirname + os.sep + "reduce.dat")

res = None
resname = None
topology = []
atoms = None
atompos = 0
nrbeads = 0
for l in open(pdb):
    if not l.startswith("ATOM"): continue
    atom = l[12:16].strip()
    cres = l[21:26]
    if cres != res:
        if atoms:
            topology.append((resname, atoms))
            nrbeads += len(ff[resname])
        res = cres
        resname = l[17:20].strip()
        assert resname in ff, l
        atoms = {}
    atoms[atom] = atompos
    atompos += 1
if atoms:
    topology.append((resname, atoms))
    nrbeads += len(ff[resname])

nratoms = atompos
assert nratoms == npy.shape[1], (nratoms, npy.shape[1])

dtype = np.dtype([
    ("coor", (np.float32, 3)),
    ("beadindex", np.uint16),
    ("beadnames", "|S4"),
])

outp = np.zeros((len(npy), nrbeads),dtype=dtype)

temp = np.zeros((len(npy), 3))
beadpos = 0
for resnr, res in enumerate(topology):
    resname, atoms = res
    ffres = ff[resname]
    for bead in ffres:
        beadindex, beadname, beadatoms, charge = bead
        temp[:] = 0
        for beadatom in beadatoms:
            assert beadatom in atoms, (resnr+1, resname, beadatom)
            atompos = atoms[beadatom]
            temp += npy[:, atompos]
        outp["coor"][:,beadpos] = temp / len(beadatoms)
        beadpos += 1

# For now, only output coordinates...
np.save(outpf, outp["coor"])
