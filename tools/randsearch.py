"""
Tool to generate random initial DOF files
syntax: python randsearch.py <bodies> <structures> [seed] [--fast] [--fix-receptor] [radius=35]
In case of two bodies:
  The bodies are placed <2*radius> Angstroms from each other, at opposite ends of the origin
  If --fix-receptor, the receptor is placed at the origin
In case of 1 or >2 bodies:
  The bodies are placed on the surface of a sphere of <radius> Angstroms
  If --fix-receptor, the first body (receptor) is placed at the origin
  The (remaining) positions are regularized, enforcing constant relative distances towards each other 
   (unless the --fast option was specified)  
"""
from __future__ import print_function

radius = 35.0

import sys
try:
  import psyco
  psyco.full()
except:
  pass  

fast = False
fix_receptor = False
try:
  i = sys.argv.index("--fast")
  sys.argv.pop(i)
  fast = True
except ValueError:
  pass  

try:
  i = sys.argv.index("--fix-receptor")
  sys.argv.pop(i)
  fix_receptor = True
except ValueError:
  pass  

try:
  i = sys.argv.index("--radius")
  radius = float(sys.argv.pop(i+1))
  sys.argv.pop(i)
except ValueError:
  pass

try:
  assert len(sys.argv) in (3, 4)
  bodies = int(sys.argv[1])
  assert bodies > 0 and bodies < 1000, bodies
  structures = int(sys.argv[2])
  assert structures > 0, structures
  seed = 1
  if len(sys.argv) > 3: seed = int(sys.argv[3])
except:
  print("syntax: python randsearch.py <bodies> <structures> [seed] [--fast] [--fix-receptor] [--radius <radius>]" , file=sys.stderr)
  raise

import random
if seed != -1:
  random.seed(seed)

rsq = radius * radius
from math import *

def adjust(x,y,z):
  w = sqrt(rsq/(x*x+y*y+z*z))
  return x*w,y*w,z*w

print("#pivot auto")
print("#centered receptor: true")
print("#centered ligands: true")
for n in range(structures):
  print("#"+str(n+1))
  p = []
  for nn in range(bodies):
    x = 2*random.random()-1 
    y = 2*random.random()-1
    z = 2*random.random()-1
    x,y,z = adjust(x,y,z)
    p.append((x,y,z))  
  if bodies == 1:
    if fix_receptor: p[0] = (0,0,0)
  elif bodies == 2 and not fix_receptor:
    x,y,z = p[0]
    p[1] = (-x,-y,-z)
  elif bodies == 3 and fix_receptor:
    p[0] = (0,0,0)
    x,y,z = p[1]
    p[2] = (-x,-y,-z)
  else:
    first = 0
    if fix_receptor: 
      p[0] = (0,0,0)
      first = 1
    while 1:  
      delta = 0
      newp = []
      for i in range(first,bodies):
	x,y,z = p[i]
	xold,yold,zold = x,y,z
	for j in range(first,bodies):
          if i==j: continue
          xx,yy,zz = p[j]
	  dx,dy,dz = xold-xx,yold-yy,zold-zz
	  dsq = dx*dx+dy*dy+dz*dz
	  dsq2 = radius/dsq
	  x += dx * dsq2
	  y += dy * dsq2
	  z += dz * dsq2
	  dx,dy,dz = x-xx,y-yy,z-zz
	x,y,z = adjust(x,y,z)
	if not fast:
          dd = x-xold,y-yold,z-zold
	  delta+=dd[0]*dd[0]+dd[1]*dd[1]+dd[2]*dd[2]
	newp.append((x,y,z))
      p[first:] = newp
      #print delta
      if fast or delta < 0.001*bodies: break
  for i in range(bodies):
    x,y,z = p[i]
    if i == 0 and fix_receptor:
      phi, ssi, rot = 0, 0, 0
    else:
      phi,ssi,rot = 2*(2*random.random()-1)*pi,2*(2*random.random()-1)*pi,2*(2*random.random()-1)*pi
    print(phi,ssi,rot,x,y,z)
