# -*- coding: utf-8 -*-
"""
Making harmonic modes for cartesian flexibility
Created on Tue Mar 26 10:15:42 2013

@author: Christina Schindler
Usage: receptor_interface.txt ligand_interface.txt receptorr.pdb ligandr.pdb
"""

import numpy as np
import sys
import re

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
        ilist = np.loadtxt(files[i],dtype=int,ndmin=1)
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
        output += '\t0\n'
    else:
        print "number of ligand hm is", len(rflex)*3
        for i, mode in enumerate(rflex):
            for j in range(3):
                output += '    '+str(3*i+j+1)+'    0.0\n'
                modes = ['0' for k in range(3*rlen)] 
                modes[3*int(mode[0]-1)+j] = '1'
                #This uncommented part is for making hydrogens flexible in conjunction with the heavy atom they are attached to
                # This was used to reduce the number of modes
                #h = find_hatoms(mode[2], mode[3])
                #if len(h) > 0:
                    #for atom in h:
                        #kk = [x[0] for x in atomlist if x[1] == mode[1] and x[2] == atom]
                        #if len(kk) == 1:
                            #modes[3*(kk[0]-1)+j] = '1'
                            
                        #else:
                            #print "ERROR", atom, "in res", mode[1], "not found"
                            #print "Check if this is consistent with your expectations"
                            #break
                            
                    #modes[3*int(mode[0]-1)+j] = '1'
                    
                #else:
                    #modes[3*int(mode[0]-1)+j] = '1'
                    
                    
                for k in range(0,len(modes),6):
                    outstring = ' '.join(modes[k:k+6])
                    output += outstring+'\n'
                
                output += '\n'
                
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
        
#-----Main--------
if parse_arg(sys.argv):   
    import subprocess
    import os
    name = sys.argv[1]
    inputfiles = sys.argv[2:]
        
    flexbeads, lenligands, atomlist = read(inputfiles)
    rcut = 3.0
    directory = os.path.split(inputfiles[1])[0]
    inp = inputfiles[0].split('rl')[0]
    while (len(flexbeads[0])+len(flexbeads[1]))*3 > 1000 and rcut > 0:
        rcut -= 0.1
        subprocess.call(['python','tools/interface.py',inputfiles[1],inputfiles[3],directory,str(rcut),'rlist-'+name,'llist-'+name])
        flexbeads, lenligands, atomlist = read(inputfiles)
        
    output = open(directory+'/flexm-'+name+'.dat','w')
    for i, beads in enumerate(flexbeads):
        tmp = gen_mode(beads, lenligands[i], name, atomlist[i])
        if i == len(flexbeads)-1:
            tmp = tmp[:-1]
            
        output.write(tmp)
        
    output.close()