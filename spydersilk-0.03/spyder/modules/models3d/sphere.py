# Copyright 2007-2009, Sjoerd de Vries
# This file is part of the Spyder module: "models3d" 
# For licensing information, see LICENSE.txt 

def spheresample(nrpoints):
  import random
  import math
  points = []
  def normpoint(p):
    size = math.sqrt(p["x"]**2+p["y"]**2+p["z"]**2)
    p["x"] /= size
    p["y"] /= size
    p["z"] /= size
  for n in range(nrpoints):
    points.append({})
    points[-1]["x"] = random.uniform(-1,1)
    points[-1]["y"] = random.uniform(-1,1)
    points[-1]["z"] = random.uniform(-1,1)
    normpoint(points[-1])
  #iters = 50 #better if you want really good positions
  iters = 10 #only this is feasable for large numbers of points (1000+)
  for i in range(iters):
    #print i
    forces = []
    for n in range(nrpoints):
      forces.append({})
      forces[-1]["x"] = 0
      forces[-1]["y"] = 0
      forces[-1]["z"] = 0
    for n in range(nrpoints):
      for nn in range(n+1,nrpoints):
        dx = points[n]["x"]-points[nn]["x"]
        dy = points[n]["y"]-points[nn]["y"]
        dz = points[n]["z"]-points[nn]["z"]
        dis = (dx**2+dy**2+dz**2)
        forcesize = 10000
        if dis > 0: forcesize = 1/dis
        forces[n]["x"] += forcesize*dx
        forces[n]["y"] += forcesize*dy
        forces[n]["z"] += forcesize*dz
        forces[nn]["x"] -= forcesize*dx
        forces[nn]["y"] -= forcesize*dy
        forces[nn]["z"] -= forcesize*dz
    for n in range(nrpoints):
      #forcesize = forces[n]["x"]**2+forces[n]["y"]**2+forces[n]["z"]**2
      #if n<10: print i, n, forcesize
      #oldx,oldy,oldz = points[n]["x"], points[n]["y"], points[n]["z"]
      points[n]["x"] += forces[n]["x"]
      points[n]["y"] += forces[n]["y"]
      points[n]["z"] += forces[n]["z"]
      normpoint(points[n])
      #newx,newy,newz = points[n]["x"], points[n]["y"], points[n]["z"]
      #if n < 10: print i,n,oldx-newx,oldy-newy,oldz-newz
      
  return points
  
