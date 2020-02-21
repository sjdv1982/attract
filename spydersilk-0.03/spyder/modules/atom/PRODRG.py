# Copyright 2008, 2009 Sjoerd de Vries
# This file is part of the Spyder module: "atom" 
# For licensing information, see LICENSE.txt 

"""
Script to contact the PRODRG server and retrieve CNS parameters
Can handle multiple small molecules, naming conflicts are automatically resolved
"""
from __future__ import print_function
import os, sys, shutil, re, math

python3 = (sys.version_info[0] == 3)
if python3:
  import http.client as httplib
  import urllib.parse as urlparse
else:
  import httplib, urlparse
import mimetypes, urllib

site = 'http://davapc1.bioch.dundee.ac.uk/'
prodrg = 'cgi-bin/prodrg/prodrg.cgi'
   
def post_multipart(pdbdata):
    """
    Contact the PRODRG server and get the CNS parameters
    """
    url = site + prodrg
    fields = (('clax_chir', ''),('clax_cha', 'WAT'), ('clax_em', 'EM'))

    urlparts = urlparse.urlsplit(url)
    host = urlparts[1]
    selector = urlparts[2]
       
    content_type = 'multipart/form-data; boundary=----abcdefgh123456789----'
    h = httplib.HTTPConnection(host)
    headers = {
      'User-Agent': 'anonymous',
        'Content-Type': content_type    
    }
    body = "coords="+pdbdata.replace(" ","+").replace("\n", "%0D%0A") 
    for f in fields:
      body = body + "&"
      body = body + f[0]
      body = body + "="
      body = body + f[1]
    h.request('POST', selector, body, headers)
    res = h.getresponse() 
    return res.read()

def submit(pdbdata):
  """
  Function to submit a PDB data with a single small molecule to the PRODRG server
  and to get CNS parameters
  The small ligand name is a three-letter 'residue' code  
  """
  
  name = ""
  for l in pdbdata.splitlines():
    if not l.startswith("ATOM"): continue
    name = l[17:20]
    break
  result = post_multipart(pdbdata)
  mode = 0
  skip = 0
  pdb, top, par = "","",""
  for l in result.split("\n"):
    ll = l.split()
    if (len(ll) > 0):
      if ll[0] == "PRODRG>" and mode == 0: mode = 1 
      if ll[0] != "PRODRG>" and mode == 1: mode = 2
      if mode == 2 and l == "<font color=yellow size=6>The PDB file (polar hydrogens)</font><br><br>":
        mode = 3
        skip = 2
      if mode > 2 and l == "</textarea>": mode = 2
      if mode == 2 and l == "<font color=yellow size=6>The CNS topology</font><br><br>":
        mode = 4
        skip = 2
      if mode == 2 and l == "<font color=yellow size=6>The CNS topology (parameters)</font><br><br>":
        mode = 5
        skip = 2
      if mode == 5 and l.upper().startswith("NBONDS"):
        mode = 6
    if skip > 0:
      skip -= 1
      continue
    if mode == 1: 
      if l.startswith("PRODRG> WARNING: atoms with same name found. Auto-renaming."):
        raise Exception("Small ligand %s contains duplicate atom names" % name)
    if mode == 3: pdb += l + "\n"
    if mode == 4: top += l + "\n"
    if mode == 5: par += l + "\n"
    if mode == 6 and l.upper().startswith("END"):
      mode = 5
  return pdb,top,par

def __matchCoor(first,second,radius=0.05):

  """
  Function to evaluate if the distance between the first and second
  coordinate is less or equal to radius. Radius has to be defined as
  radius**2 to save some CPU time. Function returns True if equal or
  smaller, else False.
  """

  d = float(0)
  for i in range(len(first)):
      d = d + (first[i]-second[i])**2

  if d <= radius: return 1
  else: return 0
  
def __cleanPRODRGoutput(prodrg_dir, tempdir, ligandname):
  
  """Remove the prodrg directory and temporary ligand pdb when done"""    
  
  shutil.rmtree(prodrg_dir)
  os.remove("%s/%s.pdb" % (tempdir,ligandname))

def __removeNBONdsFromParfile(parfile):
  
  """
  Function to remove the non-bonded parameter defenitions from the PRODRG
  parameter files. They interfere with the HADDOCK CNS defined ones and mess
  up stuff.
  """  
  
  if not parfile == None:
    
    newpar = []; read = True
    for line in parfile:
      if line.startswith('NBONds') and read:
        read = False
      elif line.startswith('END') and not read:
        read = True
      elif read:
        newpar.append(line)  
      else:
        pass
        
    return "".join(newpar)
    
  else:
    return parfile  
  
def runLocalPRODRG(pdbdata, ligandname, haddockerror, prodrgdir, tempdir):
  
  """
  Function to run a local installation of the PRODRG software package 
  on the pdbdata. Local version does not run a EM using GROMACS like the 
  PRODRG server used to do. We give return the same pdb as provided as input.
  PRODRG renames atoms, we make a mapping dictionary an rename everything back 
  in the topology and parameter files.
  """
  currdir = os.getcwd()
  os.chdir(tempdir)
  
  # Run PRODRG
  ligandname = ligandname.strip()
  tempfile = file('%s.pdb' % ligandname,'w')
  tempfile.write(pdbdata)
  tempfile.close()
  
  cmd = "%s/run_prodrg.py %s PDBELEM\n" % (prodrgdir, "%s.pdb" % ligandname)
  os.system(cmd)
  
  top, par = "",""
  prodrg_dir = "%s/%s_prodrg" % (tempdir,ligandname)
  if os.path.isdir(prodrg_dir):
    prodrg_top = "%s/%s.TOP" % (prodrg_dir, ligandname)
    prodrg_par = "%s/%s.PAR" % (prodrg_dir, ligandname)
    prodrg_fin = "%s/%s-FIN.pdb" % (prodrg_dir, ligandname)
    if os.path.isfile(prodrg_top) and os.path.isfile(prodrg_par) and os.path.isfile(prodrg_fin):
      top_file = open(prodrg_top, 'r')
      par_file = open(prodrg_par, 'r')
      fin_file = open(prodrg_fin, 'r')
      
      # Making a mapping of the atoms names in the original pdb input agains the
      # renamed atoms in the prodrg modified input pdb file.
      inresids = {}
      mapping = {}
      missing = []

      for line in pdbdata.split('\n'):
        line = line.strip()
        if line.startswith('ATOM') or line.startswith('HETATM'):
          inresids[(float(line[30:38]), float(line[38:46]), float(line[46:54]))] = line[12:16].strip() 

      for line in fin_file.readlines():
        line = line.strip()
        if line.startswith('ATOM') or line.startswith('HETATM'):
          coor = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
          found= False
          # If the coordinates match exactly than map directly
          if coor in inresids:
            mapping[line[12:16].strip()] = inresids[coor]
            found = True
          # Else evaluate if they are close together and match.
          else:
            for incoor in inresids:
              if __matchCoor(coor, incoor):
                mapping[line[12:16].strip()] = inresids[incoor]
                found = True
                break
          if not found:
            residue = line[12:16].strip()
            if not re.match('\dH.*|H.*', residue): missing.append(residue)
      
      # Final check if all heavy atoms are mapped
      if len(missing):
        __cleanPRODRGoutput(prodrg_dir, tempdir, ligandname)
        raise haddockerror("Mismatch between the atoms in the database and the supplied atoms for ligand %s:  missing: %s" 
                            % (ligandname, str(missing)))
      
      # Use the atom mapping to rename atoms in the topology files      
      proclist = list(mapping.keys())
      #print(mapping)
      for line in top_file.readlines():
        splitted = line.split()
        newline = line
        for atom in proclist:
          if atom in splitted:
            newline = re.sub(r'\s\%s\b' % atom, r' %s' %mapping[atom], newline, 1)
        top += newline
      
      # No need to remap the parameter file as it uses TYPE naming but we do need
      # to correct for NBONds parameter block by removing it all together.
      par = __removeNBONdsFromParfile(par_file.readlines())
      
      top_file.close()
      par_file.close()
      fin_file.close()
    else:
      
      # Something went wrong. Parse the PRODRG .out file for warnings
      errormsg= ""
      prodrg_out = "%s/%s_PRODRG.out" % (prodrg_dir, ligandname)
      if os.path.isfile(prodrg_out):
        out_file = open(prodrg_out, 'r')
        
        read = False
        for line in out_file.readlines():
          line = line.strip()
          if line.startswith('ERRDRG>'):
            read = False
          if read:
            errormsg += "%s\n" % line
          if line == 'PRODRG> PDB mode detected.': 
            read = True  
        
        out_file.close()
      
      __cleanPRODRGoutput(prodrg_dir, tempdir, ligandname)
      raise haddockerror("Unable to generate topology for ligand %s. PRODRG did not create the required output:\n%s" % (ligandname, errormsg))
  else:
    __cleanPRODRGoutput(prodrg_dir, tempdir, ligandname)
    raise haddockerror("Unable to generate topology for ligand %s. PRODRG did not create any output" % ligandname)
  
  __cleanPRODRGoutput(prodrg_dir, tempdir, ligandname)
  
  #Just move back to the orginal execution directory for the sake of consistency.    
  os.chdir(currdir)
  
  return pdbdata,top,par

def substitute(top,par,substituted):
  disallowed = ["H","HC","HA","C","CCIS","CH1E","CH2E","CH3E","CH2G","CH2P","C5W","CW","CR1E","C5","CRH","CR1H","CR1W","CRHH","CF","CY","CY2","N","NR","NH1","NH2","NH3","NC2","O","O1P","O3R","OC","OH1","OHP","P","SH1E","SM","S",]
  types = []
  for l in top.splitlines():
    ll = l.strip().split()
    if len(ll) == 0 or ll[0].upper() != "MASS": continue
    types.append(ll[1])
  for c0 in list(range(65,65+26))+list(range(48,58)):
    c = chr(c0)
    if c in substituted: continue
    ok = True
    for t in types:
      if len(t) <= 2: sub = t[0]+c
      else: sub = t[0]+c+t[2:]
      if sub in disallowed:
        ok = False
        break
    if ok == True: break 
  else:
    raise Exception("Too many naming conflicts between multiple PRODRG files!")
    
  substituted.add(c)
  newtop, newpar = "",""
  for l in top.splitlines():    
    ll = l.strip().split()
    if len(ll) > 0:
     typeis = l.find("TYPE=")
     if typeis > -1: typeis += len(l[typeis+5:]) - len(l[typeis+5:].strip())
     if ll[0].upper() == "MASS": 
       newtype = ll[1][0] + c
       if len(ll[1]) > 2: newtype += ll[1][2:]
       l = l[:l.index(ll[1])] + newtype + l[l.index(ll[1])+len(newtype):]
     elif typeis > -1:
       l = l[:typeis+6] + c + l[typeis+7:]
    newtop += l + '\n'
  for l in par.splitlines():
    places = []    
    if l.startswith("eval"):
      key = l[28:]
      if key.startswith("BOND"):
        places = (35,41)
      elif key.startswith("ANGLE"):
        places = (36,42,48)
      elif key.startswith("IMPR"):
        places = (35,41,47,53)
      elif key.startswith("DIHE"):
        places = (35,41,47,53)
    if l.startswith("NONBONDED"):
      places = (12,)
    for p in places:
      l = l[:p] + c + l[p+1:]    
    newpar += l + '\n'
  return newtop, newpar, substituted
    
if __name__ == '__main__':
  
  pdb = open('sah.pdb').read()
  pdbdata,top,par = runLocalPRODRG(pdb, 
                                   ligandname='ligand', 
                                   prodrgdir='/home/software/PRODRG', 
                                   tempdir='/home/enmr/services/HADDOCK/spyder/Spyder/atom/tt')
  print(top)
