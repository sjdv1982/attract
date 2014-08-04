import subprocess,sys
import argparse
import os,re, numpy
from math import *
import euler
def read_pdb(f):
  ret = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    x,y,z = (float(f) for f in (l[30:38],l[38:46],l[46:54]))
    ret.append((x,y,z))
  return ret

def get_euler(pdb1, pdb2):
  atoms1 = read_pdb(pdb1)
  atoms2 = read_pdb(pdb2)

  if len(atoms1) != len(atoms2):
    raise Exception("Different atom numbers: %s: %d, %s: %d" % (pdb1, len(atoms1), pdb2, len(atoms2)))

  atoms1 = numpy.array(atoms1)
  atoms2 = numpy.array(atoms2)

  return euler.euler(atoms1, atoms2)

  
attractdir = os.path.abspath(os.path.split(__file__)[0])
if len(attractdir) == 0: attractdir = './'
parser = argparse.ArgumentParser()
parser.add_argument('dirpath',help="path to HADDOCK structures")
parser.add_argument('receptor')
parser.add_argument('ligand')
parser.add_argument('--score',help="HADDOCK score file")
parser.add_argument('--ens',nargs=5,help='--ens receptorconf ligandconf HADDOCK_complex_file ATTRACT_ensemblefile_r ATTRACT_ensemblefile_l')
parser.add_argument('--water',action='store_const',const='w',default='',help='water refined structures')
args = parser.parse_args()

print("#pivot auto\n#centered receptor: false\n#centered ligands: false")
score = []
if args.score is not None:
  score = open(args.score).readlines()
  
ens = []
conf_r = []
conf_l = []
convert = []
if args.ens is not None:
  ens = open(args.ens[2]).readlines()
  if int(args.ens[0]) > 1:
    conf_r = open(args.ens[3]).readlines()
    
  if int(args.ens[1]) > 1:
    conf_l = open(args.ens[4]).readlines()
    
  for i in range(int(args.ens[0])):
    for j in range(int(args.ens[1])):
      convert.append((i+1,j+1))
  
for i in range(1,201):
  print("#"+str(i))
  number = str(i)+args.water
  if args.score is not None:
    filename = score[i-1].split()[0]
    filename = filename[:-5]
    number = filename.split('_')[-1]
    
  data = open(args.dirpath+'/complex_'+number+'.pdb').readlines()
  data = [line for line in data if 'ATOM' in line]
  out1 = open('tmprec.pdb','w')
  out2 = open('tmplig.pdb','w')
  for line in data:
    if ' A ' in line:
      out1.write(line)
    elif ' B ' in line:
      out2.write(line)
    else:
      print line
      
  out1.close()
  out2.close()
  receptorfile = args.receptor
  ligandfile = args.ligand
  myens = (1,1)
  if args.ens is not None:
    used = ens[int(number.replace('w',''))-1]
    used = re.sub('\D','',used)
    myens = convert[int(used)-1]
    if len(conf_r)> 0:
      receptorfile = conf_r[myens[0]-1].replace('\n','')
      
    if len(conf_l) > 0:
      ligandfile = conf_l[myens[1]-1].replace('\n','')
  
  ret = get_euler(receptorfile,'tmprec.pdb')
  ret = ["%.5f" % x if fabs(x) > 1e-5 else 0 for x in ret ]
  if len(conf_r) > 0:
    print(str(myens[0])+' '+" ".join(ret))
  else:
    print(" ".join(ret))
    
  ret = get_euler(ligandfile,'tmplig.pdb')
  ret = ["%.5f" % x if fabs(x) > 1e-5 else 0 for x in ret ]
  if len(conf_l) > 0:
    print(str(myens[1])+' '+" ".join(ret))
  else:
    print(" ".join(ret))
  
  