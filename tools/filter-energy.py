import sys
from _read_struc import read_struc
header,structures = read_struc(sys.argv[1])
energy_threshold = float(sys.argv[2])
rev = False
if len(sys.argv) > 3:
    assert len(sys.argv) == 4
    assert sys.argv[3].startswith("--rev")
    rev = True

stnr = 0
for h in header: print h
for r, struc in enumerate(structures):
    l1, l2 = struc
    for ll in l1:
        if ll.startswith("## Energy:"):
            e = float(ll[10:])
            break
    if rev:
        if e < energy_threshold: continue
    else:
        if e > energy_threshold: continue
    stnr += 1
    print "#"+str(stnr)
    print "##"+str(r) + " => filter-energy"
    try:
        for l in l1: print l
        for l in l2: print l
    except IOError:
        sys.exit()
