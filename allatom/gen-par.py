import sys
from math import sqrt

transfile = sys.argv[1]

types = {}
for line in open(transfile):
  ll = line.split()
  nr = int(ll[0])
  epsilon = float(ll[1])  
  sigma = float(ll[2])
  types[nr] = (epsilon, sigma)
  
crosstypes = {}
for nr1 in types:
  eps1, sig1 = types[nr1]
  for nr2 in types:
    eps2, sig2 = types[nr2]
    if nr1 == nr2:
      abc, rbc = 4*eps1, sig1
    else:
      abc = 4 * sqrt(eps1*eps2)
      rbc = (sig1+sig2)/2
    crosstypes[nr1,nr2] = (abc,rbc) 
    crosstypes[nr2,nr1] = (abc,rbc) 
    
abcs = [[]]
rbcs = [[]]
for i in range(1,99+1):
  cabcs = []
  crbcs = []
  for j in range(1,99+1):
    abc, rbc = 0,0
    if (i,j) in crosstypes:
      abc, rbc = crosstypes[i,j]
    cabcs.append(abc)
    crbcs.append(rbc)
  abcs.append(cabcs)
  rbcs.append(crbcs)
print 12, 99, 0, 0
for crbcs in rbcs:
  for rbc in crbcs: print "%.3f" % rbc,
  print
for cabcs in abcs:
  for abc in cabcs: print "%.3f" % abc,
  print
for i in range(1,99+1): 
  for j in range(1,99+1): 
    print 1,
  print
print




