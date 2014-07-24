# -*- coding: utf-8 -*-
"""
Making harmonic modes for cartesian flexibility
Created on Tue Mar 26 10:15:42 2013

@author: Christina Schindler
Usage: receptor_interface.txt ligand_interface.txt receptorr.pdb ligandr.pdb
"""

import numpy as np
import os,sys
import re
import interface

def parse_arg(arguments):
    if len(arguments) < 6 or not len(arguments)%2 == 0 :
        print "Incorrect number of parameters!"
        print " Usage: cartmode name receptor_interface.txt receptorr.pdb ligand_interface.txt ligandr.pdb"
        print "An arbirtrary number of ligands with corresponding interface list can be specified"
        sys.exit(1)
    else:
        return True
        
def read(files):
    # Read interface residues from files
    beadslist = []
    lengthlist = []
    atomlist = []
    for i in range(0,len(files),2):
        ilist = np.loadtxt(files[i],dtype=int)
        # Make list of corresponding beads
        data = []
        for line in open(files[i+1]):
	    tmp = line.replace('-',' -')
            list = tmp.split()
            if len(list) > 0 and list[0] == 'ATOM':
                data.append((int(list[1]),int(list[4]), list[2], list[3]))
        
        beads = []
        if len(ilist) > 0:
            for res in ilist:
                tmp =  [a for a in data if a[1] == res]# and len(re.findall('\AH', a[2])) == 0]
                beads.extend(tmp)
          
        #beads.sort()
        beadslist.append(beads)
        lengthlist.append(len(data))
        atomlist.append(data)

    return beadslist, lengthlist, atomlist  
    
def clean(mystring):
    for i in range(11,len(mystring),12):
        mystring = mystring[0:i]+'\n'+mystring[i+1:]
        
    return mystring
    
def gen_mode(rflex, rlen, name, atomlist):
    output = ''
    if len(rflex) == 0:
        output += '-1\n'
    else:
        output += '-1\n'
        for i, mode in enumerate(rflex):
	    output += str(mode[0])+'\n'

                
    return output

def find_hatoms(atomname, resname):
    #return []
    if atomname == 'N' and resname == 'PRO':
        return []
        
    elif atomname == 'N':
        return ['HN']
        
    elif atomname == 'OH':
        return ['HH']
        
    elif atomname == 'OG1':
        return ['HG1']
        
    elif atomname == 'ND2':
        return ['HD21', 'HD22']
        
    elif atomname == 'NE2' and not resname == 'HIS':
        return ['HE21', 'HE22']
      
    elif atomname == 'ND1' and resname == 'HIS':
        return ['HD1']
      
    elif atomname == 'NE2' and resname == 'HIS':
	return ['HE2']
      
    elif atomname == 'SG' and resname == 'CYS':
	return ['HG']
        
    elif atomname == 'OG':
        return ['HG']
        
    elif atomname == 'NE1':
        return ['HE1']
        
    elif atomname == 'NZ':
        return ['HZ1','HZ2','HZ3']
        
    elif atomname == 'NH1':
        return ['HH11','HH12']
        
    elif atomname == 'NH2':
        return ['HH21','HH22']
        
    elif atomname == 'NE':
        return ['HE']
        
    else:
        return []

def run(args,rcut=3.0):
  name = args[0]
  inputfiles = args[1:]
  flexbeads, lenligands, atomlist = read(inputfiles)
  directory = os.path.split(inputfiles[1])[0]
  if len(directory) == 0: directory = '.'
  inp = inputfiles[0].split('rl')[0]
  while (len(flexbeads[0])+len(flexbeads[1]))*3 > 1000 and rcut > 0:
    rcut -= 0.1
    interface.make_interface(directory+'/'+name+'.pdb', directory, name,rcut)
    flexbeads, lenligands, atomlist = read(inputfiles)
        
  output = open(directory+'/flexm-'+name+'.dat','w')
  outstring = ''
  for i, beads in enumerate(flexbeads):
    tmp = gen_mode(beads, lenligands[i], name, atomlist[i])
    if i == len(flexbeads)-1:
      if tmp == '\t0\n':
	tmp = ''
	
    outstring += tmp
            
  outstring = outstring[:-1]
  output.write(outstring)
  output.close()  
    
#-----Main--------
if __name__ == "__main__":
  if parse_arg(sys.argv):   
      run(sys.argv[1:])
        

