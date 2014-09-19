"""
This script generates the restraints necessary to preserve secondary structure in flexible docking
Advanced elastic network based on unbound protein structure
"""
import sys
import numpy as np
import math
import os
import re
import cPickle as pickle
currdir = os.path.abspath(os.path.split(__file__)[0])
allatomdir = currdir + "/../allatom"
sys.path.append(allatomdir)
from parse_cns_top import *

def make_interfacelist(ilist, pdb):
     # Read interface residues from files
    data = open(ilist).readlines()
    data = [l for l in data if not l[0]=='#']
    if len(data) == 0:
      return [],[]

    rlist = np.loadtxt(ilist,dtype=int)
    # Make list of corresponding atoms
    ratoms = []
    receptorid = []
    data = open(pdb).readlines()
    data = [x for x in data if 'ATOM' in x]
    for count,line in enumerate(data):
        tmp = line.replace('-',' -')
        list = tmp.split()
        if int(list[4]) in rlist:
	  ratoms.append(count+1)
        receptorid.append((list[3],list[4],list[2],list[1]))
                           
    return ratoms, receptorid  

def make_model(filelist,residues,presidues,nlistcut=30):
    directory, pdb, interface, name = filelist[:4]
    offset = 0
    nlist_all = pickle.load(open(pdb+'.tree',"rb"))
    if len(filelist) > 4 and not filelist[4] == '--strong':
        offset = int(filelist[4])
        
    interfacelist, atomid = make_interfacelist(interface, pdb)
    nbonds = []
    #search for neighbors and bonds
    for atom in interfacelist:
        natom = atomid[atom-1][2]
        resn = atomid[atom-1][0]
        resi = int(atomid[atom-1][1])
        a = residues[resn.lower()].copy()
        bonds = [ (b[0].upper(),b[1].upper()) for b in a.bonds]
        bonds = bonds + [('N','HT1'),('N','HT2'),('N','HT3'),('C','OXT')]
        angles = [ (ang[0].upper(),ang[2].upper()) for ang in a.get_angles()]
        angles = angles + [('HT1','HT2'),('HT2','HT3'),('HT1','HT3'),('HT1','CA'),('HT2','CA'),('HT3','CA'),
	      ('O','OXT'),('CA','OXT')]
        #Find and write bonds
        nlist = nlist_all[atom-1]
        if not nlistcut == 30:
	  nlist = nlist[:nlistcut+1]
	  
        for n,dd in nlist:
            n2 = atomid[n-1][2]
            resn2 = atomid[n-1][0]
            resi2 = int(atomid[n-1][1])
            #check = [x for x in nbonds if (x[0] == atom and x[1] == n) or (x[1] == atom and x[0] == n)]
            #if len(check):
	      #continue
	    if resi2 == resi and ((natom, n2) in bonds or (n2, natom) in bonds):              
		nbonds.append((atom, n, 1, dd)) #Write 1-2 bonds
                        
	    elif resi2 == resi and ((natom,n2) in angles or (n2,natom) in angles):
	      nbonds.append((atom, n, 2, dd)) #Write 1-3 bonds within amino acid
                  
            elif resi2 == resi + 1 or resi2 == resi - 1:
		#write babckbone bonds and angles with neighboring amino acids
		if resi2 == resi + 1 and (natom,n2) in [('C','N'),("O3'",'P')]:
		  nbonds.append((atom, n, 1, dd)) #Write 1-2 bonds
		      
		elif resi2 == resi - 1 and (natom,n2) in [('N','C'),('P',"O3'")]:
		  nbonds.append((atom, n, 1, dd)) #Write 1-2 bonds
		      
		elif resi2 == resi + 1 and (natom, n2) in [('O','N'),('C','HN'),('C','H'),('C','CA'),('CA','N'),("C3\'",'P'),("O3\'",'O1P'),("O3\'",'O2P'),("O3\'","O5\'")]:
		  nbonds.append((atom, n, 2, dd)) #Write 1-3 bonds
		      
		elif resi2 == resi - 1 and (natom, n2) in [('N','O'),('CA','C'),('H','C'),('HN','C'),('N','CA'),('O1P',"O3\'"),('O2P',"O3\'"),("O5\'","O3\'"),('P',"C3\'")]:
		  nbonds.append((atom, n, 2, dd)) #Write 1-3 bonds
		      
		else:
		  nbonds.append((atom, n, 3, dd))# Write non-bonded interaction: at least 1-4, atoms connected by bonds and angles are excluded
		      
	    elif resn == 'CYS' and resn2 == 'CYS' and not resi == resi2:
	      if natom == 'SG' and n2 == 'SG':
		# check for disulphurbridge
		if dd < 3.0:
		  nbonds.append((atom, n, 1, dd)) #Write 1-2 bonds
			
		else:
		  nbonds.append((atom, n, 3, dd)) # Write non-bonded interaction: at least 1-4, atoms connected by bonds and angles are excluded
		  
	      elif (natom,n2) in [('SG','CB'),('CB','SG')]:
		# check for disulphurbridge
		nbonds.append((atom, n, 2, dd)) #Write 1-3 bonds
			
			
	      else:
		nbonds.append((atom, n, 3, dd)) # Write non-bonded interaction: at least 1-4, atoms connected by bonds and angles are excluded
            
            else:
	      nbonds.append((atom, n, 3, dd)) # Write non-bonded interaction: at least 1-4, atoms connected by bonds and angles are excluded
            

    return nbonds, pdb, name, offset, atomid  
  
def write_output(nbonds, pdb, name, offset, atomid,strong=1.0):
    sel = []
    output = os.path.splitext(pdb)[0]+'_'+name+'.txt'
    out = open(output,'w')
    if len(nbonds) == 0:
      out.write('dummy 1 1')

    for bond in nbonds:
        #Write selection
        res1 = atomid[bond[0]-1]
        res2 = atomid[bond[1]-1]
        sel1 = res1[0]+res1[1]+'_'+res1[2]+res1[3]
        sel2 = res2[0]+res2[1]+'_'+res2[2]+res2[3]
        if not sel1 in sel:
            out.write(sel1+' 1 '+str(bond[0]+offset)+'\n')
            sel.append(sel1)
        
        if not sel2 in sel:
            out.write(sel2+' 1 '+str(bond[1]+offset)+'\n')
            sel.append(sel2)
        
    out.write('\n')
      
    check_duplicate = []
    for bond in nbonds:
        #Write restraints section
        res1 = atomid[bond[0]-1]
        res2 = atomid[bond[1]-1]
        sel1 = res1[0]+res1[1]+'_'+res1[2]+res1[3]
        sel2 = res2[0]+res2[1]+'_'+res2[2]+res2[3] 
        if (sel1,sel2) in check_duplicate or (sel2,sel1) in check_duplicate:
	  continue
        if bond[2] == 1:# and len(re.findall('\AH',res1[2])) == 0 and len(re.findall('\AH',res2[2])) == 0:
            #Write 1-2 bonds with Hinsen model elastic bonds
            out.write(sel1+' '+sel2+' 4 '+str(bond[3])+' '+str(1000.0*strong)+'\n')
        
        elif bond[2] == 2:# and len(re.findall('\AH',res1[2])) == 0 and len(re.findall('\AH',res2[2])) == 0:
            #Write 1-3 bonds to fix angles
            out.write(sel1+' '+sel2+' 4 '+str(bond[3])+' '+str(100.0*strong)+'\n')
            
        #Write restraints for other parts of the elastic network with force constant
        #and minimum distance restraints (type 3)
        elif bond[2] == 3:
            sigma1 = 3.5
            sigma2 = 3.5
            #ignore hydrogens in steric repulsion
            if 'H' in res1[2] and len(re.findall('[CNO]H', res1[2])) == 0:
                sigma1 = 0.0
                
            if 'H' in res2[2] and len(re.findall('[CNO]H', res2[2])) == 0:
                sigma2 = 0.0
                
            mindist = math.sqrt(sigma1*sigma2)
            if mindist > 0: 
                if bond[3] < mindist:
                    mindist = bond[3]
                
                out.write(sel1+' '+sel2+' 5 '+str(mindist)+' '+str(1.0*strong)+'\n')  
        
        check_duplicate.append((sel1,sel2))
        
    out.close()    

def make_restraints(topology, *args):
    residues, presidues = topology
    nbonds, pdb, name, offset, atomid = make_model(args,residues,presidues)
    cut = 31
    while len(nbonds) > 9000 and cut > 0:
      cut -= 5
      nbonds, pdb, name, offset, atomid = make_model(args,residues,presidues,cut)
      
    if '--strong' in args:
      write_output(nbonds,pdb,name,offset,atomid,10.0)
    else:
      write_output(nbonds,pdb,name,offset,atomid)
  
#Main
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        ilist, pdb = sys.argv[1:]
        atomlist, pdb2, pdbid = make_interfacelist(ilist,pdb)
        outfile = ''
        if 'A' in pdb:
            outfile = os.path.dirname(pdb)+'/ilistA.txt'
            
        elif 'B' in pdb:
             outfile = os.path.dirname(pdb)+'/ilistB.txt'           
            
        out = open(outfile,'w')
        for atom in atomlist:        
            out.write(str(atom)+'\n')
            
        out.close()
        
    else: 
        topology = parse_stream(open(sys.argv[1]))
        make_restraints(topology,*sys.argv[2:])
