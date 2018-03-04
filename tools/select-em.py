import sys
from _read_struc import read_struc

def read_enestruc(datfile):
    header,structures = read_struc(datfile)
    structures = list(structures)
    energies = []
    for l1,l2 in structures:
      for ll in l1:
        if ll.startswith("## Energy:"):
          ee = ll[10:].strip()
          if ee.startswith("nan"):
            e = 99999999999999
          else:
            e = float(ee)
          energies.append(e)
    energies = energies[:4]
    energy = sum(energies)/len(energies)
    return header, energy, structures

runs = int(sys.argv[1])
iter_refinements = int(sys.argv[2])
data = []
header = None
for n in range(runs):
    filename = "dock-refine%d-it%d-topstruc.dat" % (n+1, iter_refinements)
    header0, energy, structures = read_enestruc(filename)
    if n == 0:
        header = header0
    else:
        assert header == header0, (n+1)
    data.append((energy, structures))
data.sort(key=lambda v:v[0], reverse=True)

for h in header: print(h)
stnr = 0
for n in range(4):
    for run in range(runs):
        energy, structures = data[run]
        s = structures[n]
        stnr += 1
        l1,l2 = s
        print("#"+str(stnr))
        print("## Run %d structure %d" % (run+1, n+1))
        print("## Energy: %.6f" % energy )
        for l in l2: print(l)

for run in range(runs):
    energy, structures = data[run]
    for n in range(4,len(structures)):
        s = structures[n]
        stnr += 1
        l1,l2 = s
        print("#"+str(stnr))
        print("## Run %d structure %d" % (run+1, n+1))
        print("## Energy: %.6f" % energy )
        for l in l2: print(l)
