# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:05:30 2013

@author: christina

iATTRACT refinement with OPLS force field
Supports ensemble docking (tested) and multi-body docking (untested)
"""

import sys, random, os, time, itertools
from math import *
from multiprocessing import Pool
from copy import deepcopy
currdir = os.path.abspath(os.path.split(__file__)[0])
if len(currdir) == 0: currdir = '.'
attractdir = currdir + "/../bin"
tooldir = currdir + "/../tools"
allatomdir = currdir + "/../allatom"
sys.path.append(allatomdir)
sys.path.append(tooldir)
sys.path.append(attractdir)
import imodes, get_restraints
import parse_cns_top
import collectlibpy as collectlib
import neighbortree

topstream = [(allatomdir + "/topallhdg5.3.pro", "oplsx"),
             (allatomdir + "/dna-rna-allatom.top", "dna-rna")
            ]

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

def read_pdb(f):
  ret1 = []
  for l in open(f):
    if not l.startswith("ATOM"): continue
    ret1.append(l)
  return ret1

#prepare input for on the fly flexible interface refinement
def prepare_input(topology,start,ligands,current,name,coor,ligandrange,ligandatoms,ensemblefiles=[],modefile=None,otf=True,noflex=[],icut=3.0):
  current = str(current)
  directorypath = os.path.split(ligands[0])[0]
  if len(directorypath) == 0: directorypath = '.'
  currligands = ligands
  if len(ensemblefiles):
    data = open(start).readlines()
    data = [line for line in data if not '#' in line]
    for nr, filename in ensemblefiles:
      currensnr = int(data[int(nr)-1].split()[0])
      ensnames = open(filename).readlines()
      currligands[int(nr)-1] = ensnames[currensnr-1].replace('\n','')

  if otf and os.path.exists('flexm-'+current+name+'.dat'):
    check = True
    restraints = []
    for ligand in currligands:
      if not os.path.exists(os.path.splitext(ligand)[0]+'_'+current+name+'.rest'):
	check = False
      else:
	restraints.append(os.path.splitext(ligand)[0]+'_'+current+name+'.rest')

    if check:
      return ('flexm-'+current+name+'.dat',restraints)

  if otf:
    imodes.make(ligands,ligandatoms,current+name,ligandrange,coor,icut,ensemblefiles,modefile,imodefile,noflex)

  else:
    restraints = []
    for i, ligand in enumerate(currligands):
      restraints.append(os.path.splitext(ligand)[0]+'_'+name+'.rest')

    return ('flexm-'+name+'.dat',restraints)


  if len(noflex) == len(ligands):
    return ("",[])

  offset = 0
  restraints = []
  for i,ligand in enumerate(currligands):
    if len(ilist):
      get_restraints.make_restraints(topology,directorypath,ligand,ilist[i],current+name,str(offset))
    else:
      get_restraints.make_restraints(topology,directorypath,ligand,current+name+'-ilist'+str(i+1)+'.txt',current+name,str(offset))

    restraints.append(os.path.splitext(ligand)[0]+'_'+current+name+'.rest')
    offset += len(read_pdb(ligand))

  return ('flexm-'+current+name+'.dat',restraints)

def prepare_input2(topology,ilist, ligands, name, args):
    directorypath = os.path.split(ligands[0])[0]
    if len(directorypath) == 0: directorypath = '.'
    ensemblefiles = []
    for i, item in enumerate(args):
      if item == '--ens':
	ensemblefiles.append((args[i+1],args[i+2]))

    imodes.make_defined(ilist,ligands,name)
    offset = 0
    for i in range(len(ligands)):
      tmp = [item[1] for item in ensemblefiles if int(item[0]) == i+1]
      if len(tmp):
	data = open(tmp[0]).readlines()
	for line in data:
	  filename = line.replace('\n','')
	  get_restraints.make_restraints(topology,directorypath,filename,ilist[i],name,str(offset))

      else:
	get_restraints.make_restraints(topology,directorypath,ligands[i],ilist[i],name,str(offset))

      offset += len(read_pdb(ligands[i]))

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
    if len(ligands) > 2:
      newargs = []
      for i, item in enumerate(args):
	if not item in ligands:
	  newargs.append(item)

      args = newargs[:1]+['partners-aa.pdb']+newargs[1:]

    for i, item in enumerate(args):
      if item == '--ens':
	ensemblefiles.append((args[i+1],args[i+2]))

      elif item == '--modes':
	modefile = args[i+1]

    imodefile, restfiles = prepare_input(topology,start,ligands,current,name,coor,ligandrange,ligandatoms,ensemblefiles,modefile,otf,noflex,icut)

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
  exists = False
  while 1:
    anr += 1

    if anr > len(sys.argv)-1: break
    arg = sys.argv[anr]
    if arg == "--torque":
      torque = "-torque"
      sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
      anr -= 1
      continue

    if arg == "--exists":
      exists= True
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

    if arg == "--top" or arg == "--topfile":
      assert os.path.exists(nextarg), nextarg
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

  for f in topfiles:
    topstream.append((f, f))
  for f, nam in topstream:
    parse_cns_top.parse_stream(open(f), nam)
  topology = parse_cns_top.residues, parse_cns_top.presidues
  totstruc = get_struc(strucfile)
  if jobsize is not None:
    chunks = ceil(totstruc/float(jobsize))
  if not chunks: sys.exit()

  args = sys.argv[2:]
  if torque == '-infinite':
    args = args + ['--rcut','1']


  scoremode = "--score" in args
  ligandrange,coor,ligandatoms = [], [], []
  if not '--vmax' in args:
    args = args + ['--vmax','2500']

  ligands = [item for item in args if '.pdb' in item]
  ligandatoms = []
  ensfiles = []
  modefile = None
  imodefile = None
  for i, item in enumerate(args):
    if item == '--ens':
      ensfiles.append((args[i+1],args[i+2]))

    elif item == '--modes':
      modefile = args[i+1]

    elif item == '--imodes':
      imodefile = args[i+1]

  for i,u in enumerate(ligands):
    ligandatoms.append(read_pdb(u))
    has_ens = None
    for nr, ensfile in ensfiles:
      if str(i+1) == nr:
	has_ens = ensfile

    if has_ens is not None:
      ensligand = open(has_ens).readlines()
      for ligandname in ensligand:
	neighbortree.make(ligandname.replace('\n',''))

    else:
      neighbortree.make(u)

  if len(ligands) > 2:
    os.system('echo '+str(len(ligands))+' > partners-aa.pdb')
    for l in ligands:
      os.system('cat '+l+' >> partners-aa.pdb')
      os.system('echo TER >> partners-aa.pdb')

  #check if interface lists have been provided and generate global imodes and restraints files from them
  if '--ilist' in args:
    otf = False
    k = args.index('--ilist')
    ligands = [item for item in args if '.pdb' in item]
    if '--' in args[k+1:k+len(ligands)+1]:
      raise ValueError("You must specify one interfacelist for each ligand --ilist ilist-1.txt ilist-2.txt []")

    check = False
    for i in range(len(ligands)):
      if not os.path.exists(args[k+1+i]):
	check = True

    if check:
      raise ValueError("At least one of the interface list files does not exist, did you generate them correctly?")

    ilist = args[k+1:k+len(ligands)+1]
    args = args[:k]+args[k+len(ligands)+1:]
    prepare_input2(topology,ilist,ligands,name,args)

  if otf:
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
    while 1:
      if exists: break
      result = collectlib.collect_next()
      if result: break
      nstruc += 1
      tmpcoor = collectlib.collect_all_coor()
      coor.append(deepcopy(tmpcoor))

  else:
    coor = [[0] for i in range(len(ligands))]
  if exists:
    for i in range(int(chunks)):
      coor.append([])

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
    com = "python2 %s/split.py %s %s %d" % (tooldir, strucfile, pat, chunks); run(com)
    runs = [(i,attract,pat,pat2,args,name,ligandrange,ligandatoms,coor[i-1],otf,noflex,interface_cutoff) for i in range(1,int(chunks)+1)]
    p = Pool(np)
    try:
      #for crun in runs: run_docking(crun) #serial execution
      p.map_async(run_docking,runs).get(99999) #parallel execution

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
    #check all files if they contain correct structures
    for i in range(1,int(chunks)+1):
      data = open(pat2+'-'+str(i)).readlines()
      if not len(data) or not len(data[-1]) > 6:
	com = "cp %s-%d %s-%d" % (pat,i,pat2,i)
	run(com)

    com = "python2 %s/join.py %s %s >> %s" % (tooldir, pat2, score, output) 
    run(com)
  finally:
    com = "rm %s-*" % pat; run(com)
    com = "rm %s-*" % pat2; run(com)
