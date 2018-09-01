from __future__ import print_function
import parmed, sys

def usage():
    print("%s\nusage: python trans_par_amber.py <parameter file> <atom type start>\n%s" % (60*"*",60*"*"), file=sys.stderr)

try:
    fil = sys.argv[1]
    start = int(sys.argv[2])
except:
    usage()
    raise

p = parmed.amber.AmberParm(fil)

dict_sigeps = {}

for a in p.atoms:
    if (a.sigma, a.epsilon) not in dict_sigeps.keys():
        dict_sigeps[(a.sigma, a.epsilon)] = set()
    dict_sigeps[(a.sigma, a.epsilon)].add(a.type)

i = start
for (s, e) in dict_sigeps:
    print('%i %.4f %.4f '%(i, e, s), end=''),
    types = dict_sigeps[(s, e)]
    for t in types:
        print(t.upper(), sep=' ', end='')
    print('')
    i+=1
