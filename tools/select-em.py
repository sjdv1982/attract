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
    energies4 = energies[:4]
    energy = sum(energies4)/len(energies4)
    return header, energy, energies, structures

runs = int(sys.argv[1])
iter_refinements = int(sys.argv[2])
data = []
header = None
for n in range(runs):
    filename = "dock-refine%d-it%d-topstruc.dat" % (n+1, iter_refinements)
    header0, energy, energies, structures = read_enestruc(filename)
    if n == 0:
        header = header0
    else:
        assert header == header0, (n+1)
    data.append((energy, energies, structures))
data.sort(key=lambda v:v[0], reverse=True)

for h in header: print(h)
stnr = 0
for n in range(4):
    for run in range(runs):
        energy, energies, structures = data[run]
        s = structures[n]
        e = energies[n]
        stnr += 1
        l1,l2 = s
        print("#"+str(stnr))
        print("## Run %d structure %d" % (run+1, n+1))
        print("## Energy: %.6f" % e )
        for l in l2: print(l)

for run in range(runs):
    energy, energies, structures = data[run]
    for n in range(4,len(structures)):
        s = structures[n]
        e = energies[n]
        stnr += 1
        l1,l2 = s
        print("#"+str(stnr))
        print("## Run %d structure %d" % (run+1, n+1))
        print("## Energy: %.6f" % e )
        for l in l2: print(l)
