# -*- coding: utf-8 -*-
"""
Find all interface residues
Created on Wed Jun 12 13:14:50 2013

@author: christina
"""
def read_file(file1):
    atomlist = []
    for line in open(file1):
        x,y,z = (float(f) for f in (line[30:38],line[38:46],line[46:54]))
        if l.startswith('ATOM'):
            atomlist.append((int(line[4:11]),int(line[22:26]),x,y,z))
            
    return atomlist

def read_struc(file1):
    atomlista = []
    atomlistb = []
    count = 0
    data = open(file1).readlines()
    data = [ x for x in data if x.startswith('ATOM')]
    for count, l in enumerate(data):
        x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
        item = (int(line[4:11]),int(line[22:26]), x, y, z)
        if count < item[0]:            
            atomlista.append(item)            
        elif count > item[0]:
            atomlistb.append(item)
          
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

def collect_contacts(a1,a2,Y,rcut):
  ilist1, ilist2 = [],[]
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
	  
  return ilist1, ilist2
	  
import math  
import numpy as np
from scipy.spatial.distance import cdist
def get_interface(a1, a2, output,output2,rcut):
    ilist1 = []
    ilist2 = []
    crd1 = [[x,y,z] for name, res, x,y,z in a1]
    crd1 = np.matrix(crd1)
    crd2 = [[x,y,z] for name, res, x,y,z in a2]
    crd2 = np.matrix(crd2)
    Y = cdist(crd1,crd2,'euclidean')
    while (len(ilist1) < 8 or len(ilist2) < 8) and countilist(ilist1,a1,ilist2,a2) < 333:
	ilist1, ilist2 = collect_contacts(a1,a2,Y,rcut)             
        rcut += 0.5
     
    while countilist(ilist1,a1,ilist2,a2) > 333:
      rcut -= 0.5
      ilist1, ilist2 = collect_contacts(a1,a2,Y,rcut)
      
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

def get_contacts(structure,rcut=3.0):
    a1, a2 = read_struc(structure)
    contacts = []
    ncontact = 0
    crd1 = [[x,y,z] for name, res, x,y,z in a1]
    crd1 = np.matrix(crd1)
    crd2 = [[x,y,z] for name, res, x,y,z in a2]
    crd2 = np.matrix(crd2)
    Y = cdist(crd1,crd2,'euclidean')
    for i in range(len(a1)):
      res1 = a1[i][1]
      for j in range(len(a2)):
	res2 = a2[j][1]
	dist = Y[i][j]
	if dist < rcut:
	  if not (res1,res2) in contacts:
	    print res1, res2
	    contacts.append((res1,res2))
	    ncontact += 1
	    
    return contacts
	  
def get_interface2(a1, a2, output,ouput2,rcut):
    ilist = []
    while len(ilist) < 8:
        rcut2 = rcut*rcut
        for atom1 in a1:
            for atom2 in a2:
                if not atom1[1] in ilist:
                    dist = (atom1[2]-atom2[2])**2+(atom1[3]-atom2[3])**2+(atom1[4]-atom2[4])**2
                    if dist < rcut2:
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
        
def check_contacts(struc,ilist1old, ilist2old,directory,name):
    a1, a2 = read_struc(struc)
    ilist1, ilist2 = [], []
    crd1 = [[x,y,z] for name1, res, x,y,z in a1]
    crd1 = np.matrix(crd1)
    crd2 = [[x,y,z] for name1, res, x,y,z in a2]
    crd2 = np.matrix(crd2)
    Y = cdist(crd1,crd2,'euclidean')
    remove = open(directory+'/'+name+'rlist-remove.txt','w')
    ilist1new, ilist2new = collect_contacts(a1,a2,Y,5.0)
    for res in ilist1old:
      if not res in ilist1new:
	remove.write(str(res)+'\n')
      else:
	ilist1.append(int(res))
	
    remove.close()
    remove2 = open(directory+'/'+name+'llist-remove.txt','w')
    for res in ilist2old:
      if not res in ilist2new:
	remove2.write(str(res)+'\n')
	
      else:
	ilist2.append(int(res))
	
    ilist1backup = ilist1[:]
    ilist2backup = ilist2[:]
    remove2.close()
    rcut = 3.0
    while (len(ilist1) < 8 or len(ilist2) < 8) and countilist(ilist1,a1,ilist2,a2) < 333:
      ilist1new, ilist2new = collect_contacts(a1,a2,Y,rcut)
      for res in ilist1new:
	if not res in ilist1:
	  ilist1.append(res)
	  
      for res in ilist2new:
	if not res in ilist2:
	  ilist2.append(res)
	  
      rcut += 0.5
     
    while countilist(ilist1,a1,ilist2,a2) > 333:
      rcut -= 0.5
      ilist1new, ilist2new = collect_contacts(a1,a2,Y,rcut)
      for res in ilist1:
	if not res in ilist1new and not res in ilist1backup:
	  ilist1.pop(ilist1.index(res))
	  
      for res in ilist2:
	if not res in ilist2new and not res in ilist2backup:
	  ilist2.pop(ilist2.index(res))
	  
	  
    output = directory+'/'+name+'rlist.txt'
    output2 = directory+'/'+name+'llist.txt'
    out = open(output, 'w')
    if len(ilist1) == 0:
        out.write('# no flexible residues')
        
    else:
        for i in ilist1:
            out.write(str(i)+'\n')
            
    out.close()      
    out2 = open(output2, 'w')
    if len(ilist2) == 0:
        out2.write('# no flexible residues')
        
    else:
        for i in ilist2:
            out2.write(str(i)+'\n')
            
    out2.close()
    
    
import os    
def make_interface(struc,directory,name,rcut=3.0):
    a1, a2 = read_struc(struc)
    get_interface(a1, a2, directory+'/'+name+'rlist.txt',directory+'/'+name+'llist.txt',rcut)
    
#Main
if __name__ == "__main__":
    import sys
    if '--ncontact' in sys.argv:
      get_contacts(sys.argv[1])
      
    elif '--continue' in sys.argv:
      struc, fileilist1, fileilist2,directory,name = sys.argv[1:6]
      ilist1 = list(np.loadtxt(fileilist1,dtype='int'))
      ilist2 = list(np.loadtxt(fileilist2,dtype='int'))
      check_contacts(struc,ilist1,ilist2,directory,name)
      
    elif len(sys.argv) == 3:
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
