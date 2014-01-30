# -*- coding: utf-8 -*-
"""
Find all interface residues
Created on Wed Jun 12 13:14:50 2013

@author: christina
"""
def read_file(file1):
    atomlist = []
    for line in open(file1):
        tmp = line.replace('-',' -')
        list = tmp.split()
        if len(list) > 0 and list[0] == 'ATOM':
            atomlist.append((int(list[1]),int(list[4]),float(list[5]),float(list[6]),float(list[7])))
            
    return atomlist

def read_struc(file1):
    atomlista = []
    atomlistb = []
    count = 0
    data = open(file1).readlines()
    data = [ x for x in data if 'ATOM' in x]
    for count, line in enumerate(data):
        tmp = line.replace('-',' -')
        l = line
        list = tmp.split()
        if count < int(list[1]):
            x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
            atomlista.append((int(list[1]),int(list[4]),x,y,z))
            
        elif count > int(list[1]):
            x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
            atomlistb.append((int(list[1]),int(list[4]),x,y,z))
          
    return atomlista, atomlistb

def countilist(ilist1, a1, ilist2, a2):
  count = 0
  for atom in a1:
    if atom[1] in ilist1:
      count += 1
      
  for atom in a2:
    if atom[1] in ilist2:
      count += 1
      
  return count

import math  
import numpy as np
from scipy.spatial.distance import cdist
def get_interface(a1, a2, output,output2,rcut):
    ilist1 = []
    ilist2 = []
    contactlist = []
    distdata = {}
    crd1 = [[x,y,z] for name, res, x,y,z in a1]
    crd1 = np.matrix(crd1)
    crd2 = [[x,y,z] for name, res, x,y,z in a2]
    crd2 = np.matrix(crd2)
    Y = cdist(crd1,crd2,'euclidean')
    while (len(ilist1) < 8 or len(ilist2) < 8) and countilist(ilist1,a1,ilist2,a2) < 333:
	for i in range(len(a1)):
	  res1 = a1[i][1]
	  for j in range(len(a2)):
	    res2 = a2[j][1]
	    dist = Y[i][j]
	    if dist < rcut:
	      if not res1 in ilist1:
		ilist1.append(res1)
                        
	      if not res2 in ilist2:
		ilist2.append(res2)
                        
        rcut += 0.5
     
    count = 0
    while countilist(ilist1,a1,ilist2,a2) > 333:
      if count%2 == 0 or count%3 == 0:
	ilist1.pop()
      else:
	ilist2.pop()
      
      count += 1
      
    out = open(output, 'w')
    if len(ilist1) == 0:
        out.write('# no flexible residues')
        
    else:
        ilist1.sort()
        for i in ilist1:
            out.write(str(i)+'\n')
            
    out.close()      
    out2 = open(output2, 'w')
    if len(ilist2) == 0:
        out2.write('# no flexible residues')
        
    else:
        ilist2.sort()
        for i in ilist2:
            out2.write(str(i)+'\n')
            
    out2.close()
        
def get_interface2(a1, a2, output,ouput2,rcut):
    ilist = []
    while len(ilist) < 8:
        rcut2 = rcut*rcut
        for atom1 in a1:
            for atom2 in a2:
                if not atom1[1] in ilist:
                    dist = (atom1[2]-atom2[2])**2+(atom1[3]-atom2[3])**2+(atom1[4]-atom2[4])**2
                    if dist < rcut2:
                        #print dist, atom1, atom2
                        ilist.append(atom1[1])
                        break
                    
                else:
                    break
                        
        rcut += 0.5
                    
    out = open(output, 'w')
    if len(ilist) == 0:
        out.write('# no flexible residues')
        out.close()
        
    else:
        for i in ilist:
            out.write(str(i)+'\n')
            
        out.close()
    
def make_interfacelist(file1, file2, directory,rcut=3.0,name1='rlist',name2='llist'):
    a1 = read_file(file1)
    a2 = read_file(file2)
    get_interface(a1, a2, directory+'/'+name1+'.txt',directory+'/'+name2+'.txt',rcut)
    #if not name1 == 'rlist':
        #import restraints_from_topology as resttop
        #import numpy as np
        #ilist = np.loadtxt(directory+'/'+name1+'.txt',dtype=int)
        #atomlist, pdb2, pdbid = resttop.make_interfacelist(ilist,file1)
        #outfile = directory+'/'+name1+'-aa.txt'
        #out = open(outfile,'w')
        #for atom in atomlist:        
            #out.write(str(atom)+'\n')
            
        #out.close()
        #ilist = np.loadtxt(directory+'/'+name2+'.txt',dtype=int)
        #atomlist, pdb2, pdbid = resttop.make_interfacelist(ilist,file2)
        #outfile = directory+'/'+name2+'-aa.txt'
        #out = open(outfile,'w')
        #for atom in atomlist:        
            #out.write(str(atom)+'\n')
            
        #out.close()
        

import os    
def make_interface(struc,directory,name,rcut=3.0):
    a1, a2 = read_struc(struc)
    get_interface(a1, a2, directory+'/'+name+'rlist.txt',directory+'/'+name+'llist.txt',rcut*rcut)
    
#Main
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        struc, directory = sys.argv[1:]
        make_interface(struc, directory)
        
    elif len(sys.argv) == 4:
        struc, directory,name = sys.argv[1:]
        make_interface(struc, directory, name)
    
    elif len(sys.argv) == 5:
        struc, directory,name,rcut = sys.argv[1:]
        make_interface(struc, directory, name,float(rcut))    
    
    elif len(sys.argv) == 7:
        file1, file2, directory, rcut, name1, name2 = sys.argv[1:]
        make_interfacelist(file1, file2, directory,float(rcut), name1, name2)
        
    else:
        file1, file2, directory, rcut = sys.argv[1:]
        make_interfacelist(file1, file2, directory,float(rcut))
