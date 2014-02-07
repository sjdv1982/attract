import subprocess
import numpy as np
import sys
import os
import glob
import random

structures = sys.argv[1]
directory = os.path.split(structures)[0]
nstruc = sys.argv[2]
lenA = int(sys.argv[3])
filterlist = []
ATTRACT='/home/ichauvot/attract'
pattern = random.randint(0,9999)
subprocess.call('python '+ATTRACT+'/tools/split.py '+structures+' tmp'+str(pattern)+' '+nstruc,shell=True)
aspfile=open('SASaspfile','w')
phefile=open('SASphefile','w') 
for i in range(int(nstruc)):
  if len(filterlist) == 10000:
    break
  
  filename = 'tmp'+str(pattern)+'-'+str(i+1)
  subprocess.call(ATTRACT+'/bin/collect '+filename+' '+directory+'protein.pdb '+directory+'peptide/peptide.pdb --ens 2 '+directory+'peptide/peptides.list > tmp'+str(pattern)+'.pdb',shell=True)
  subprocess.call(['rm',filename])
  p = subprocess.Popen('/home/ichauvot/attract/bin/shrake tmp'+str(pattern)+'.pdb 1.4', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  asp = 0
  phe = 0
  for line in p.stdout.readlines():
      if ' '+str(lenA+1)+' ' in line:
       tmp = line.split()
       asp = float(tmp[-1])
      
      elif ' '+str(lenA+4)+' ' in line:
       tmp = line.split()
       phe = float(tmp[-1])
	
  retval = p.wait()
  if phe < 75.0 and asp > 100.0:
#  if phe < 1000.0 and asp > 1.0:
    filterlist.append(i+1)
  aspfile.write('%i\n'%asp)
  phefile.write('%i\n'%phe)
from math import *
from _read_struc import read_struc  

header,structures = read_struc(sys.argv[1])
structures = list(structures)

for h in header: print h

stnr = 0
st2nr = 0
for s in structures:
  stnr += 1
  if not stnr in filterlist:
    continue
  else:
    st2nr += 1
    print "#"+str(st2nr)
    print "## "+str(stnr) + " => filter"
    l1,l2 = s
    for l in l1: print l
    for l in l2: print l
