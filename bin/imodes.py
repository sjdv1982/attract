""" -------------------------
Get interface for iATTRACT calculation and make imodes file
supports ensemble und multiple ligands

Author: CEM Schindler
31/7/2014
-----------------------------"""

import math  
import numpy as np
from scipy.spatial.distance import cdist

def validate(ilist,unbounds):
  check = False
  count = 0
  for i in range(len(ilist)):
    if len(ilist[i]) < 8:
      check = True
      
    flex = [l for l in unbounds[i] if int(l.split()[4]) in ilist[i]]
    count += len(flex)
   
  if count > 333 and check == True:
    check = False
    
  return check

def countflex(ilist,unbounds):
  check = True
  count = 0
  for i in range(len(ilist)):
    flex = [l for l in unbounds[i] if int(l.split()[4]) in ilist[i]]
    count += len(flex)
   
  if count < 333:
    check = False
    
  return check  

def collect_contacts(unbounds,ligandrange,Y,rcut):
  ilist = []
  allunbounds = []
  for b in unbounds:
    allunbounds.extend(b)
    ilist.append([])
   
  access = []
  for i in range(len(unbounds)):
    for j in range(i+1,len(unbounds)):
      for ii in range(ligandrange[i][0],ligandrange[i][1]):
	for jj in range(ligandrange[j][0],ligandrange[j][1]):
	  access.append((ii,jj))
    
  contacts = [ item for item in access if Y[item[0]][item[1]] < rcut] 
  for ii,jj in contacts:
    for i in range(len(unbounds)):
      if ii >= ligandrange[i][0] and ii < ligandrange[i][1]:
	line = allunbounds[ii]
	res1 = int(line.split()[4]) 
	if not res1 in ilist[i]:
	  ilist[i].append(res1)
	  
      if jj >= ligandrange[i][0] and jj < ligandrange[i][1]:
	line2 = allunbounds[jj]
        res2 = int(line2.split()[4])
	if not res2 in ilist[i]:
	  ilist[i].append(res2)

  for i in range(len(ilist)):
    ilist[i].sort()
    
  return ilist
	  
  
def get_interface(ligandcoor,ligandpdbs,ligandrange,rcut):
    ilist = [[],[]]
    ligandcoor = np.matrix(ligandcoor)
    Y = cdist(ligandcoor,ligandcoor,'euclidean')
    while validate(ilist,ligandpdbs):
	ilist = collect_contacts(ligandpdbs,ligandrange,Y,rcut)
        rcut += 0.5
     
    while countflex(ilist,ligandpdbs):
      rcut -= 0.5
      ilist = collect_contacts(ligandpdbs,ligandrange,Y,rcut) 
      
    return ilist

def read_pdb(f):
  ret1 = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    ret1.append(l)
  return ret1

def make(ligands,ligandatoms, name,ligandrange,coor,thresh=3.0,ensfiles=[],modefile=None,imodefile=None):
  interface = get_interface(coor, ligandatoms,ligandrange,thresh)
  imodestring = ''
  for i in range(len(interface)):
    imodestring += '-1\n'
    data = [j+1 for j, line in enumerate(ligandatoms[i]) if int(line.split()[4]) in interface[i]]
    for item in data:
      imodestring += str(item)+'\n'
      
    out = open(name+'-ilist'+str(i+1)+'.txt','w')
    if len(interface[i]) == 0:
      out.write('#no flexible residues')
    for item in interface[i]:
      out.write(str(item)+'\n')

    out.close()
    
  imodestring = imodestring[:-1]
  newimode = open('flexm-'+name+'.dat','w')  
  newimode.write(imodestring)
  newimode.close()

def make_defined(ifilelist,ligands,name):
  interface = []
  for fn in ifilelist:
    #TODO check if this is correct
    interface.append(np.loadtxt(fn,dtype=int))

  ligandatoms = []
  for i in ligands:
    ligandatoms.append(read_pdb(i))
    
  for i in range(len(interface)):
    imodestring += '-1\n'
    data = [j+1 for j, line in enumerate(ligandatoms[i]) if int(line.split()[4]) in interface[i]]
    for item in data:
      imodestring += str(item)+'\n'
      
  imodestring = imodestring[:-1]
  newimode = open('flexm-'+name+'.dat','w')  
  newimode.write(imodestring)
  newimode.close()
  
#Main
if __name__ == "__main__":
    import sys
    make(sys.argv[1],sys.argv[2:-1],sys.argv[-1])
