from __future__ import print_function
import sys

def usage():
    print("%s\nusage: python trans_par.py <parameter file> <atom type start>\n%s" % (60*"*",60*"*"), file=sys.stderr)

try:
    fil = open(sys.argv[1])
    start = int(sys.argv[2])
except:
    usage()
    raise

types = []

for l in fil:
    l = l.lstrip().upper().rstrip()
    if not l.startswith("NONBONDED"): continue
    ll = l.split()
    v = ll[2:4]
    for t in types:
        if t[1] == v:
            t[0].append(ll[1])
            break
    else:
        types.append([[ll[1]], v])

for tnr,t in enumerate(types):
    print(start+tnr, end= " ")
    print(t[1][0],t[1][1], end= " ")
    for typ in t[0]:
        print(typ, end= " ")
    print()
