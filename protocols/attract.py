import sys, random, os, time
from math import *
import subprocess, threading

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
  #checking if the process has ended doesn't work for slow file systems
  if scoremode:
    fstruc = get_energy(f)
  else:
    fstruc = get_struc(f)
  return fstruc == nstruc

processes = []
def run(command, thread=False):
  print command  
  def func():
    p = subprocess.Popen([command], shell=True)
    processes.append(p)
    p.wait()    
  if thread:
    t = threading.Thread(target=func)
    t.start()
  else:
    func()  

np = 1
output = None
anr = 0
existing = None
jobsize = None
chunks = None
attractversion = ""
while 1:
  anr += 1

  if anr > len(sys.argv)-1: break  
  arg = sys.argv[anr]
  if arg == "--infinite":
    attractversion = "-infinite"
    sys.argv = sys.argv[:anr] + sys.argv[anr+1:]
    anr -= 1
    continue  
  
  if anr >= len(sys.argv)-1: break
  arg, nextarg = sys.argv[anr],sys.argv[anr+1]

  if arg.startswith("--exist"):
    existing = nextarg
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

queue = [None for n in range(np)]
nstrucs = [None for n in range(np)]

attractdir0 = os.path.split(sys.argv[0])[0]
tooldir = attractdir0 + "/../tools"
attractdir = attractdir0 + "/../bin"

attract = attractdir + "/attract" + attractversion
strucfile = sys.argv[1]

if jobsize is None and chunks is None:
  raise ValueError("You must specify --jobsize <value> or --chunks <value>")
if jobsize is not None and chunks is not None:
  raise ValueError("You must specify --jobsize <value> OR --chunks <value>, not both!")

totstruc = get_struc(strucfile)
if jobsize is not None:
  chunks = ceil(totstruc/float(jobsize))
if not chunks: sys.exit()

args = sys.argv[2:]
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
  done = 0
  current = 1
  processes = []
  while 1:
    for vnr in range(np):
      v = queue[vnr]
      if v is None: continue
      if finished(v, nstrucs[vnr]): 
	done += 1
	if done == chunks: break
	queue[vnr] = None
    if done == chunks: break

    free = [n for n,v in enumerate(queue) if v is None]
    if len(free) == 0 or current == chunks+1:
      time.sleep(2)
      continue

    q = free[0]
    inp = "%s-%d" % (pat, current)
    outp = "%s-%d" % (pat2, current)
    queue[q] = outp
    nstrucs[q] = get_struc(inp)
    com = " ".join([attract,inp]+args) + " > %s" % outp
    if existing is not None:
      ef = "%s-%d" % (existing, current)
      if os.path.exists(ef):
        com = "cp %s %s" % (ef, outp)
    run(com, thread=True)
    current += 1

  o = open(output, "w")
  print >> o, "## Command line arguments: " + " ".join([attract,strucfile]+args)
  o.close()
  score = ""
  if scoremode:
    score = "--score"  
  com = "python %s/join.py %s %s >> %s" % (tooldir, pat2, score, output) 
  run(com)  
finally:
  for p in processes:
    if p.returncode is None: #process is still running
      print "KILL"
      p.kill()
  com = "rm %s-*" % pat; run(com)
  com = "rm %s-*" % pat2; run(com)
