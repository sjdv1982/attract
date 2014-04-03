# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:05:30 2013

@author: christina
"""

import sys, random, os, time
from math import *
from multiprocessing import Pool

def get_energy(f):
  if not os.path.exists(f): return 0
  ret = 0
  f0 = open(f)
  for l in f0.readlines():
    if l.lstrip().startswith("Energy:"): ret += 1
  return ret
    
def get_struc(f):  
  if not os.path.exists(f): return 0
  ret = 0
  f0 = open(f)
  for l in f0.readlines()[-100:]:
    if not l.startswith("#"): continue
    try:
      v = int(l[1:])
    except ValueError:
      continue
    if v > ret: ret = v  
  f0.close()
  return ret

def finished(f, nstruc):
  if scoremode:
    fstruc = get_energy(f)
  else:
    fstruc = get_struc(f)
  return fstruc == nstruc

def run(command):
  print command
  os.system(command)  
 
def read_file(file1):
    a1 = []
    for line in open(file1):
        list = line.split()
        if len(list) > 0 and list[0] == 'ATOM':
            a1.append(list[1])  
            
    return len(a1)

import subprocess
#prepare input for on the fly flexible interface refinement
def prepare_input(start,pdbA,pdbB,current,name,attracttools):
  current = str(current)
  directorypath = os.path.split(pdbA)[0]
  if len(directorypath) == 0: directorypath = '.'
  subprocess.call(attracttools+'/../bin/collect '+start+' '+pdbA+' '+pdbB+' > '+directorypath+'/'+current+name+'.pdb',shell=True)
  subprocess.call(['python',attracttools+'/interface.py',directorypath+'/'+current+name+'.pdb',directorypath,current+name])
  count = 0
  if not os.path.exists(directorypath+'/'+current+name+'rlist.txt'):
    out = open(directorypath+'/'+current+name+'rlist.txt','w')
    out.write('#no flexible residues')
    out.close()
    count += 1
  
  if not os.path.exists(directorypath+'/'+current+name+'llist.txt'):
    out = open(directorypath+'/'+current+name+'llist.txt','w')
    out.write('#no flexible residues')
    out.close()
    count += 1
    
  if count == 2:
    return ("","","")
  
  subprocess.call(['python',attracttools+'/cartmode.py',current+name,directorypath+'/'+current+name+'rlist.txt',pdbA,directorypath+'/'+current+name+'llist.txt',pdbB])
  subprocess.call(['python',attracttools+'/get_restraints.py',directorypath,pdbA,directorypath+'/'+current+name+'rlist.txt',current+name])
  lenA = read_file(pdbA)
  subprocess.call(['python',attracttools+'/get_restraints.py',directorypath,pdbB,directorypath+'/'+current+name+'llist.txt',current+name,str(lenA)])
  subprocess.call(['rm',directorypath+'/'+current+name+'.pdb'])
  return (directorypath+'/flexm-'+current+name+'.dat',os.path.splitext(pdbA)[0]+'_'+current+name+'.txt',os.path.splitext(pdbB)[0]+'_'+current+name+'.txt')

    
#prepare input for run with global interface list
def prepare_input2(pdbA,pdbB,rlist,llist,name,attracttools):
  directorypath = os.path.split(pdbA)[0]  
  if len(directorypath) == 0: directorypath = '.'
  subprocess.call(['python',attracttools+'/cartmode.py',name,rlist,pdbA,llist,pdbB])
  subprocess.call(['python',attracttools+'/get_restraints.py',directorypath,pdbA,rlist,name])
  lenA = read_file(pdbA)
  subprocess.call(['python',attracttools+'/get_restraints.py',directorypath,pdbB,llist,name,str(lenA)])
    
def run_docking(datain):

    current,attract,pat,pat2,args,name,otf = datain
    start = pat+'-'+str(current)
    outp = pat2+'-'+str(current)
    imodefile = ''
    restfile1 = ''
    restfile2 = ''
    directory = os.path.split(args[1])[0]
    if len(directory) == 0: directorypath = '.'
	attracttools = os.path.split(attract)[0]+'/../tools'
    if otf:
	imodefile, restfile1, restfile2 = prepare_input(start,args[1],args[2],current,name,attracttools)
    
    else:
	imodefile = directory+'/flexm-'+name+'.dat'
	restfile1 = os.path.splitext(args[1])[0]+'_'+name+'.txt'
	restfile2 = os.path.splitext(args[2])[0]+'_'+name+'.txt'
    
    if imodefile == '':
	com = "cp %s %s" %(start,outp)
	print "No flexible residues for %s, skipping..." % start
	run(com)
	return
    
    com = " ".join([attract,start]+args+['--imodes',imodefile,'--rest',restfile1,'--rest',restfile2]) + " > %s " % outp
    run(com)


np = 1
output = None
anr = 0
torque = ""
jobsize = None
chunks = None
name = None
otf = True
while 1:
  anr += 1

  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]
  if arg == "--torque":
    torque = "-torque"
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    anr -= 1
    continue
  
  if arg == "--infinite":
    torque = "-infinite"
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    anr -= 1
    continue
  
  if anr >= len(sys.argv)-1: break
  arg, nextarg = sys.argv[anr],sys.argv[anr+1]
  
  if arg == "--name":
    name = nextarg
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 1
    continue
  
  if arg == "-np" or arg == "--np":
    try:
      np = int(nextarg)
      if np <= 0: raise ValueError
    except ValueError: 
      raise ValueError("Invalid number of processors: %s" % nextarg)
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 1
    continue
  
  if arg == "-jobsize" or arg == "--jobsize":
    try:
      jobsize = int(nextarg)
      if jobsize <= 0: raise ValueError
    except ValueError: 
      raise ValueError("Invalid jobsize: %s" % nextarg)
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 1
    continue
  
  if arg == "-chunks" or arg == "--chunks":
    try:
      chunks = int(nextarg)
      if chunks <= 0: raise ValueError
    except ValueError: 
      raise ValueError("Invalid chunks: %s" % nextarg)
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 1
    continue
  
  if arg == "--output":
    output = nextarg
    sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
    anr -= 1

if output is None:
  raise ValueError("You must specify an output file with --output")

if name is None:
  raise ValueError("You must specify a name for the run with --name")

attractdir0 = os.path.split(sys.argv[0])[0]
if len(attractdir0) == 0: attractdir0 = '.'
tooldir = attractdir0 + "/../tools"
attractdir = attractdir0 + "/../bin"

attract = attractdir + "/attract" + torque
strucfile = sys.argv[1]
if otf and not jobsize == 1:
  jobsize = 1
if jobsize is None and chunks is None:
  raise ValueError("You must specify --jobsize <value> or --chunks <value>")
if jobsize is not None and chunks is not None:
  raise ValueError("You must specify --jobsize <value> OR --chunks <value>, not both!")

totstruc = get_struc(strucfile)
print totstruc
if jobsize is not None:
  chunks = ceil(totstruc/float(jobsize))
if not chunks: sys.exit()

args = sys.argv[2:]
if torque == '-infinite':
  args = args + ['--rcut','1']
#check if interface lists have been provided and generate global imodes and restraints files from them
if '--ilist' in args:
  otf = False
  k = args.index('--ilist')
  if '--' in args[k+1] or '--' in args[k+2]:
    raise ValueError("You must specify one interfacelist for each ligand --ilist rlist.txt llist.txt")
  
  if not (os.path.exists(args[k+1]) and os.path.exists(args[k+2])):
    raise ValueError("At least one of the interface list files does not exist, did you generate them correctly?")
  
  rlist = args[k+1]
  llist = args[k+2]
  prepare_input2(args[1],args[2],rlist,llist,name,tooldir)
  args = args[:k]+args[k+3:]
  
scoremode = "--score" in args
while 1:
  pat = "tmp%d" % random.randint(1,99999)
  pat2 = "tmp%d" % random.randint(1,99999)
  if (pat == pat2): continue
  if os.path.exists("%s-1" % pat): continue
  if os.path.exists("%s-1" % pat2): continue
  break  

try:
  com = "python %s/split.py %s %s %d" % (tooldir, strucfile, pat, chunks); run(com)
  runs = [(i,attract,pat,pat2,args,name,otf) for i in range(1,int(chunks)+1)]
  p = Pool(np)
  try:
      p.map_async(run_docking,runs).get(99999)
      p.close()
  except KeyboardInterrupt:
      p.terminate()
      print "You cancelled the program!"
      sys.exit(1)
  finally:
      p.join()
      
  o = open(output, "w")
  o.write("## Command line arguments: " + " ".join([attract,strucfile]+args))
  o.write('\n')
  o.close()
  score = ""
  if scoremode:
    score = "--score"  
  com = "python %s/join.py %s %s >> %s" % (tooldir, pat2, score, output) 
  run(com)
finally:
  com = "rm %s-*" % pat; run(com)
  com = "rm %s-*" % pat2; run(com)
