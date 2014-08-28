import numpy
import sys

pdb1 = sys.argv[1]
pdb2 = sys.argv[2]

def read_pdb(pdb):
  ret = {}
  for l in open(pdb).readlines():
    if l[:4] != "ATOM" or l[13:16] != "CA ": continue
    resid = l[22:26]
    assert resid not in ret, resid
    x,y,z = float(l[30:38]),float(l[38:46]),float(l[46:54])
    ret[resid] = numpy.array((x,y,z))
  return ret
    
cA = read_pdb(pdb1)
cB = read_pdb(pdb2)
assert len(cB.keys()) == len(cB.keys())
for k in cA: assert k in cB

vec = numpy.array((0.0,0.0,0.0))
first = True
for k1 in cA:
  a1A = cA[k1]
  a1B = cB[k1]
  v1 = a1A - a1B
  for k2 in cA:
    if k2 is k1: continue
    a2A = cA[k2]
    a2B = cB[k2]
    v2 = a2A - a2B
    cross = numpy.cross(v1,v2)
    if not first:    
      dp = numpy.dot(vec,cross)
      assert(abs(dp) > 0.95)
      if dp < 0: cross *= -1    
    #print "%.3f %.3f %.3f" % tuple(cross/numpy.linalg.norm(cross))
    vec += cross
    first = False
vec /= numpy.linalg.norm(vec) 
print "%.6f %.6f %.6f" % tuple(vec)
    
    
  
  
  
