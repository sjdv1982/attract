# -*- coding: utf-8 -*-
"""
Prepare allatom structures for docking
Created on Thu Nov  7 15:32:17 2013

@author: christina
"""

import sys
import os
import subprocess

def read_file(filename):
  data = []
  pdblines = open(filename).readlines()
  for i,line in enumerate(pdblines):
    list = line.split()
    if len(list) > 0 and list[0] == 'ATOM':
      data.append((list[2],int(list[4]),i))
      
  return data, pdblines

def check(unbound, bound,sort=False):
  atomlistu, pdbu = read_file(unbound)
  atomlistb, pdbb = read_file(bound)
  if not len(atomlistu) == len(atomlistb):
    print "ERROR in length of pdb files"
    sys.exit(1)
    
  for i in range(len(atomlistu)):
    atom1 = atomlistu[i][0]
    res1 = atomlistu[i][1]
    atom2 = atomlistb[i][0]
    if not atom1 == atom2:
      if sort:
	for j in range(i+1,i+10):
	  atom3 = atomlistu[j][0]
	  res3 = atomlistu[j][1]
	  if res3 == res1 and atom3 == atom2:
	    atomlistu[i], atomlistu[j] = atomlistu[j], atomlistu[i]
	    pdbu[i], pdbu[j] = pdbu[j], pdbu[i]
	    break
	  
	else:
	  print "ATOM not found", atom1, atom2
	  print atom1, atom2
	  print unbound, bound
	  sys.exit(1)
	    
      else:
	print "ERROR: different atoms detected at line",i
	print atom1, atom2
	print unbound, bound
	sys.exit(1)
	
    if ("HE2" in pdbu[i] or "HD1" in pdbu[i]) and 'XXX' in pdbu[i]:
      if (not "HE2" in pdbb[i]) and (not "HD1" in pdbb[i]):
	print pdbb[i]
	sys.exit(1)
	
      else:
	print "Check HIS", unbound
	tmp1 = pdbu[i-1]
	tmp1 = tmp1.replace("-"," -")
	tmp2 = pdbb[i-1]
	tmp2 = tmp2.replace("-"," -")
	tmp3 = pdbb[i]
	tmp3 = tmp3.replace("-"," -")
	x1 = float(tmp1.split()[5])
	y1 = float(tmp1.split()[6])
	z1 = float(tmp1.split()[7])
	x2 = float(tmp2.split()[5])
	y2 = float(tmp2.split()[6])
	z2 = float(tmp2.split()[7])
	x3 = float(tmp3.split()[5])
	y3 = float(tmp3.split()[6])
	z3 = float(tmp3.split()[7])
	dx = x3-x2
	dy = y3-y2
	dz = z3-z2
	x = "%.3f"%(x1 + dx)
	y = "%.3f"%(y1 + dy)
	z = "%.3f"%(z1 + dz)
	prex,postx = x.split('.')
	xout = ''
	for k in range(8-len(prex)):
	  xout += ' '
    
	xout += x
	yout = ''
	prey, posty = y.split('.')
	for k in range(7-len(postx)-len(prey)):
	  yout += ' '
	  
	yout += y
	zout = ''
	prez, postz = z.split('.')
	for k in range(7-len(posty)-len(prez)):
	  zout += ' '
	  
	zout += z
	tmp2 = pdbu[i][:26]+xout+yout+zout+pdbu[i][54:]
	print pdbu[i]
	print pdbu[i-1]
	print tmp2
	pdbu[i] = tmp2
    
    if "HG" in pdbu[i] and "CYS" in pdbu[i] and 'XXX' in pdbu[i]:
      if not 'XX' in pdbb[i] and "HG" in pdbb[i]:
	print "Check CYS", unbound
	tmp1 = pdbu[i-1]
	tmp1 = tmp1.replace("-"," -")
	tmp2 = pdbb[i-1]
	tmp2 = tmp2.replace("-"," -")
	tmp3 = pdbb[i]
	tmp3 = tmp3.replace("-"," -")
	x1 = float(tmp1.split()[5])
	y1 = float(tmp1.split()[6])
	z1 = float(tmp1.split()[7])
	x2 = float(tmp2.split()[5])
	y2 = float(tmp2.split()[6])
	z2 = float(tmp2.split()[7])
	x3 = float(tmp3.split()[5])
	y3 = float(tmp3.split()[6])
	z3 = float(tmp3.split()[7])
	dx = x3-x2
	dy = y3-y2
	dz = z3-z2
	x = "%.3f"%(x1 + dx)
	y = "%.3f"%(y1 + dy)
	z = "%.3f"%(z1 + dz)
	prex,postx = x.split('.')
	xout = ''
	for k in range(8-len(prex)):
	  xout += ' '
	  
	xout += x
	yout = ''
	prey, posty = y.split('.')
	for k in range(7-len(postx)-len(prey)):
	  yout += ' '
	  
	yout += y
	zout = ''
	prez, postz = z.split('.')
	for k in range(7-len(posty)-len(prez)):
	  zout += ' '
	  
	zout += z
	tmp2 = pdbu[i][:26]+xout+yout+zout+pdbu[i][54:]
	print pdbu[i]
	print pdbu[i-1]
	print tmp2
	pdbu[i] = tmp2
	
    if "HG" in pdbu[i] and "CYS" in pdbu[i] and 'XXX' in pdbb[i] and "HG" in pdbb[i]:
      tmp2 = pdbu[i][:27]+' XXX '+pdbu[i][27:]
      pdbu[i] = tmp2
	
	
  subprocess.call(['cp',unbound,unbound+'.save'])
  out = open(unbound,'w')
  for line in pdbu:
    out.write(line)
    
  out.close()
  

#path to structural data
pathname = sys.argv[1]
dirlist = [ name for name in os.listdir(pathname) if os.path.isdir(os.path.join(pathname, name)) ]
old_dir = os.getcwd()
os.chdir(pathname)
ligands = ['A','B']
calc = True
for directory in dirlist:
  print directory
  if directory == '1DQJ':
    calc = True
    
  for ligand in ligands:
    unbound = directory+'/'+directory+ligand+'-unbound.pdb'
    unboundaa =directory+'/'+directory+ligand+'-unbound-aa.pdb' 
    bound = directory+'/'+directory+ligand+'-refe.pdb'
    boundaa = directory+'/'+directory+ligand+'-refe-aa.pdb'
   # if os.path.exists(unboundaa) and os.path.exists(boundaa):
   #   continue
    if os.path.exists(unbound) and os.path.exists(bound) and calc:
      os.chdir(old_dir+'/allatom')
      subprocess.call(['python','pqreduce-notermini-refepatches.py',old_dir+'/'+pathname+'/'+unbound,'oplsx.trans','topallhdg5.3.pro',old_dir+'/'+pathname+'/'+bound])
      os.chdir(old_dir)
      #subprocess.call(['python','tools/clean_pdb.py',pathname+'/'+unboundaa])
      #subprocess.call(['python','tools/clean_pdb.py',pathname+'/'+boundaa])
      os.chdir(old_dir+'/'+pathname)
      check(unboundaa,boundaa,True)
      subprocess.call('grep -v XXX '+unboundaa+' > tmp.pdb',shell=True)
      subprocess.call('grep -v XXX '+boundaa+' > tmpref.pdb',shell=True)
      os.chdir(old_dir)
      subprocess.call('python pdbtools/pdb_atom_renumber.py '+pathname+'/tmp.pdb > '+pathname+'/'+unboundaa,shell=True)
      subprocess.call('python pdbtools/pdb_atom_renumber.py '+pathname+'/tmpref.pdb > '+pathname+'/'+boundaa,shell=True)
      check(pathname+'/'+unboundaa,pathname+'/'+boundaa)
      subprocess.call(['cp',pathname+'/'+unboundaa,pathname+'/'+directory+'/'+directory+ligand+'.pdb'])
      os.chdir(pathname)
      
    else:
      print "Files not found in directory",directory
      print unbound, bound
      break