import sys, argparse
import numpy as np
import collectlibpy as collectlib
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
import rmsdlib
from scipy.spatial.distance import cdist

########################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('rec_u_npy')
parser.add_argument('lig_u_npy')
parser.add_argument('rec_b')
parser.add_argument('lig_b')
parser.add_argument('--cutoff', default=5 ,type=float)

args = parser.parse_args()
########################

cutoffsq = args.cutoff*args.cutoff

#unboundfiles = []
#
#for n in range(3, len(sys.argv), 2):
#  unboundfiles.append(sys.argv[n])
#  boundfiles.append(sys.argv[n+1])
boundfiles = [args.rec_b, args.lig_b]
bounds = [rmsdlib.read_pdb(f) for f in boundfiles]

unboundfiles = [args.rec_u_npy, args.lig_u_npy]
unbounds = [np.load(f) for f in unboundfiles]

rec = unbounds[0]
if rec.ndim == 2:
    unbounds[0] = rec[None,:,:]
    rec = unbounds[0]
else:
    assert len(rec) == len(unbounds[1]), "different nb of structures in rec and lig"
unboundsizes = [a.shape[1] for a in unbounds]

boundsizes = [len(list(p.atoms())) for p in bounds]
for boundsize, boundfile, unboundsize, unboundfile in zip(
        boundsizes, boundfiles, unboundsizes, unboundfiles
  ):
    assert boundsize == unboundsize, (boundfile, unboundfile, boundsize, unboundsize)

bound_hmask = rmsdlib.build_hydrogenmask(bounds)
for p in bounds:
  p.remove_hydrogens()

contacts = rmsdlib.build_contactfilters(bounds, args.cutoff)

f1 = sys.stdout
lig = unbounds[1]
nstruc = len(lig)
Nat = sum(unboundsizes)
coor0 = np.zeros((Nat,3), dtype=rec.dtype)
coor0[:unboundsizes[0]] = rec[0]

for struc in range(nstruc):
    coor = coor0
    if len(rec) != 1:
        coor[:unboundsizes[0]] = rec[nstruc]
    coor[unboundsizes[0]:] = lig[struc]
    coor = np.compress(bound_hmask, coor, axis=0)
    fcount = 0
    #TODO: make faster by pdist
    for filter1, filter2 in contacts:
        coor1 = np.take(coor, filter1, axis=0)
        coor2 = np.take(coor, filter2, axis=0)
        distances = cdist(coor1, coor2,'sqeuclidean')
        mindistance = distances.min()
        if mindistance < cutoffsq:
          fcount += 1
    f1.write("%.2f" % (float(fcount)/len(contacts)) +'\n')
