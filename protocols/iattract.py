# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:05:30 2013

@author: christina

iATTRACT refinement with OPLS force field
TODO make multi body and ensemble docking iattract refinement
"""

import sys, random, os, time
from math import *
from multiprocessing import Pool
currdir = os.path.abspath(os.path.split(__file__)[0])
if len(currdir) == 0: currdir = '.'
attractdir = currdir + "/../bin"
tooldir = currdir + "/../tools"
allatomdir = currdir + "/../allatom"
sys.path.append(allatomdir)
sys.path.append(tooldir)
sys.path.append(attractdir)
import imodes, get_restraints
import collectlibpy as collectlib
topology = allatomdir+'/topallhdg5.3.pro'
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
  
def read_pdb(f):
  ret1 = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    ret1.append(l)
  return ret1

#prepare input for on the fly flexible interface refinement
def prepare_input(start,ligands,current,name,coor,ligandrange,ligandatoms,ensemblefiles=[],modefile=None,otf=True,noflex=[],icut=3.0):
  current = str(current)
  directorypath = os.path.split(ligands[0])[0]
  if len(directorypath) == 0: directorypath = '.'
  currligands = ligands
  if len(ensemblefiles):
    data = open(start).readlines()
    # TODO check if this is the right format for dat files
    data = [line for line in data if not '#' in line]
    for nr, filename in ensemblefiles:
      currensnr = int(data[int(nr)-1].split()[0])
      ensnames = open(filename).readlines()
      currligands[int(nr)-1] = ensnames[currensnr-1]
      
  if otf and os.path.exists(directorypath+'/flexm-'+current+name+'.dat'):
    check = True
    restraints = []
    for ligand in currligands:
      if not os.path.exists(os.path.splitext(ligand)[0]+'_'+current+name+'.txt'):
	check = False
      else:
	restraints.append(os.path.splitext(ligand)[0]+'_'+current+name+'.txt')
    
    if check:
      return (directorypath+'/flexm-'+current+name+'.dat',restraints)
    
  if otf:
    imodes.make(ligands,ligandatoms,current+name,ligandrange,coor,icut,ensemblefiles,modefile,imodefile)
    
  else:
    restraints = []
    for i, ligand in enumerate(currligands):
      restraints.append(os.path.splitext(ligand)[0]+'_'+name+'.txt')
      
    return (directorypath+'/flexm-'+name+'.dat',restraints)
  
  count = 0
  if len(noflex):
    for ligand in noflex:
      out = open(directorypath+'/'+current+name+'-ilist'+str(ligand)+'.txt','w')
      out.write('#no flexible residues')
      out.close()	
	
  for i,ligand in enumerate(ligands):
    if not os.path.exists(directorypath+'/'+current+name+'-ilist'+str(i+1)+'.txt'):
      out = open(directorypath+'/'+current+name+'-ilist'+str(i+1)+'.txt','w')
      out.write('#no flexible residues')
      out.close()
      count += 1
  
    
  if count == len(ligands):
    return ("",[])
  
  offset = 0
  restraints = []
  for i,ligand in enumerate(currligands):
    if len(ilist):
      get_restraints.make_restraints([topology,directorypath,ligand,ilist[i],current+name,str(offset)])
    else:
      get_restraints.make_restraints([topology,directorypath,ligand,directorypath+'/'+current+name+'-ilist'+str(i+1)+'.txt',current+name,str(offset)])
    
    restraints.append(os.path.splitext(ligand)[0]+'_'+current+name+'.txt')
    offset += read_file(ligand)

  return (directorypath+'/flexm-'+current+name+'.dat',restraints)

def prepare_input2(ilist, ligands, name, args):
    directorypath = os.path.split(ligands[0])[0]
    if len(directorypath) == 0: directorypath = '.'
    ensemblefiles = []
    for i, item in enumerate(args):
      if item == '--ens':
	ensemblefiles.append((args[i+1],args[i+2]))
    
    imodes.make_defined(ilist,ligands,name)
    offset = 0
    for i in range(len(ligands)):
      tmp = [item[1] for item in ensfiles if int(item[0]) == i+1]
      if len(tmp):
	data = open(tmp[0]).readlines()
	for line in data:
	  filename = line.replace('\n','')
	  get_restraints.make_restraints([topology,directorypath,filename,directorypath+'/'+name+'-ilist'+str(i+1)+'.txt',name,str(offset)])
	  
      else:
	get_restraints.make_restraints([topology,directorypath,ligands[i],directorypath+'/'+name+'-ilist'+str(i+1)+'.txt',name,str(offset)])
	
      offset += read_file(ligands[i])
    
def run_docking(datain):
    current,attract,pat,pat2,args,name,ligandrange,ligandatoms,coor,otf,noflex,icut = datain
    start = pat+'-'+str(current)
    outp = pat2+'-'+str(current)
    imodefile = ''
    restfile1 = ''
    restfile2 = ''
    directory = os.path.split(args[1])[0]
    if len(directory) == 0: directory = '.'
    ensemblefiles = []
    modefile = None
    ligands = [item for item in args if '.pdb' in item]
    for i, item in enumerate(args):
      if item == '--ens':
	ensemblefiles.append((args[i+1],args[i+2]))
	
      elif item == '--modes':
	modefile = args[i+1]
	
    imodefile, restfiles = prepare_input(start,ligands,current,name,coor,ligandrange,ligandatoms,ensemblefiles,modefile,otf,noflex,icut)
    
    if imodefile == '':
	com = "cp %s %s" %(start,outp)
	print "No flexible residues for %s, skipping..." % start
	run(com)
	return
    
    rest = []
    for r in restfiles: rest += ['--rest',r]
    com = " ".join([attract,start]+args+['--imodes',imodefile]+rest) + " > %s " % outp
    run(com)

if __name__ == "__main__":
  np = 1
  output = None
  anr = 0
  torque = ""
  jobsize = 1
  chunks = None
  name = None
  ilist = []
  noflex = []
  interface_cutoff = 3.0
  otf = True
  topfiles = []
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

    if arg == "--icut":
      interface_cutoff = float(nextarg)
      sys.argv = sys.argv[:anr] + sys.argv[anr+2:]
      anr -= 1
      continue
  
    if arg == "--top":
      topfiles.append(nextarg)
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
  
    if arg == "--noflex":
      try:
	noflex.append(int(nextarg))
	if noflex[-1] <= 0 or noflex[-1] > 2: raise ValueError
      except ValueError: 
	raise ValueError("Invalid selection for ligand without flexibility: %s" % nextarg)
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


  attract = attractdir + "/attract" + torque
  strucfile = sys.argv[1]
  if jobsize is None and chunks is None:
    raise ValueError("You must specify --jobsize <value> or --chunks <value>")
  if jobsize is not None and chunks is not None:
    raise ValueError("You must specify --jobsize <value> OR --chunks <value>, not both!")

  if len(topfiles):
    import subprocess
    subprocess.call(['rm','/tmp/topology.top'])
    for topo in topfiles:
      subprocess.call('cat '+topo+' >> /tmp/topology.top',shell=True)
    
    topology = '/tmp/topology.top'
    
  totstruc = get_struc(strucfile)
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
    ligands = [item for item in args if '.pdb' in item]
    if '--' in args[k+1:k+len(ligands)]:
      raise ValueError("You must specify one interfacelist for each ligand --ilist ilist-1.txt ilist-2.txt []")
  
    check = False
    for i in range(len(ligands)):
      if not os.path.exists(args[k+i]):
	check = True
    
    if check:
      raise ValueError("At least one of the interface list files does not exist, did you generate them correctly?")
  
    ilist = args[k+1:k+len(ligands)]
    args = args[:k]+args[k+len(ligands):]
    prepare_input2(ilist,ligands,name,args)
  
  scoremode = "--score" in args
  ligandrange,coor,ligandatoms = [], [], []
  if otf:
    ligands = [item for item in args if '.pdb' in item]
    ligandatoms = []
    ensfiles = []
    modefile = None
    imodefile = None
    for u in ligands:
      ligandatoms.append(read_pdb(u))
   
    for i, item in enumerate(args):
      if item == '--ens':
	ensfiles.append((args[i+1],args[i+2]))

      elif item == '--modes':
	modefile = args[i+1]

      elif item == '--imodes':
	imodefile = args[i+1]
      
    initargs = [strucfile]+ligands
    if modefile: initargs += ["--modes", modefile]
    if imodefile: initargs += ["--imodes", imodefile]
    for nr, ensfile in ensfiles:
      initargs += ["--ens", nr, ensfile]
    
    collectlib.collect_init(initargs)
    ligandrange = []
    start0 = 0
    for i in collectlib.ieins[:len(ligands)]:
      ligandrange.append((start0,i))
      start0 = i
  
    nstruc = 0
    coor = []
    while 1:
      result = collectlib.collect_next()
      if result: break
      nstruc += 1
      coor.append(collectlib.collect_all_coor())

  if otf and not len(coor) == int(chunks):
    raise ValueError("For on the fly interface docking we need to determine the interface of all structures separately, please specify --jobsize 1")
    
  while 1:
    pat = "tmp%d" % random.randint(1,99999)
    pat2 = "tmp%d" % random.randint(1,99999)
    if (pat == pat2): continue
    if os.path.exists("%s-1" % pat): continue
    if os.path.exists("%s-1" % pat2): continue
    break  

  try:
    com = "python %s/split.py %s %s %d" % (tooldir, strucfile, pat, chunks); run(com)
    runs = [(i,attract,pat,pat2,args,name,ligandrange,ligandatoms,coor[i-1],otf,noflex,interface_cutoff) for i in range(1,int(chunks)+1)]
    p = Pool(np)
    try:
      p.map_async(run_docking,runs).get(99999)
      
    except KeyboardInterrupt:
      p.terminate()
      print "You cancelled the program!"
      sys.exit(1)
      
    finally:
      p.close()
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
