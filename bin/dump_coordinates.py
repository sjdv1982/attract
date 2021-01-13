#collects coordinates of selected ligand atoms and dumps them to a .npy file
# syntax: dump_coordinates <.dat file> <ligand PDB file> <output file> <nr_atoms> [atom indices] [arguments to collect]
# nr_atoms can be -1 (followed by no atom indices), in which case all atoms are collected

import sys, os
datfile=sys.argv[1]
ligpdbfile=sys.argv[2]
outfile=sys.argv[3]
if len(sys.argv) > 4:
  nr_dump_atoms = int(sys.argv[4])
else:
  nr_dump_atoms = -1
if nr_dump_atoms == -1:
  dump_atoms = None
  rest = sys.argv[5:]
else:
  dump_atoms = [int(v) for v in sys.argv[5:5+nr_dump_atoms]]
  rest = sys.argv[5+nr_dump_atoms:]

import numpy
import collectlibpy as collectlib

initargs = [datfile, ligpdbfile, ligpdbfile] + rest
collectlib.collect_init(initargs)
result = collectlib.collect_next()
coor = collectlib.collect_coor_raw()
assert len(coor) % 6 == 0, len(coor)
coorsize = len(coor) / 2
natom = len(coor) / 6
if dump_atoms is None:
  nr_dump_atoms = natom

maxstruc = 10**5
allcoor = numpy.zeros(shape=(maxstruc,3*nr_dump_atoms))

mask = []
for n in range(coorsize):
  mask.append(False)
for n in range(natom):
  sel = False
  if dump_atoms is None or n+1 in dump_atoms: sel=True
  for n in range(3): mask.append(sel)

mask = numpy.array(mask, dtype=bool)
allcoor[0,:] = coor[mask]
nstruc = 1

while 1:
  result = collectlib.collect_next()
  if result: break
  coor = collectlib.collect_coor_raw()
  allcoor[nstruc,:] = coor[mask]
  nstruc += 1
  if nstruc == maxstruc:
    maxstruc = int(maxstruc * 1.2 + 0.99999) #grow allcoor with 20%
    allcoor.resize(maxstruc,3*nr_dump_atoms)
allcoor = allcoor.reshape((len(allcoor), -1, 3))
numpy.save(outfile, allcoor[:nstruc])