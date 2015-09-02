"""
Converts CNS/HADDOCK .tbl file into ATTRACT restraint file
Author: Sjoerd de Vries, TU Munich

Syntax: tbl2attract <.tbl file> {[--pdb <PDB file>]} {[--mapping <mapping files>] [options]
 
 PDBs must be reduced, and for every PDB file, a mapping file must be supplied
 The first PDB will be "segid A" in the .tbl file, the second PDB "segid B", etc.
"""

import sys

import argparse
parser = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("tblfile",help="CNS/HADDOCK .tbl file to convert")
parser.add_argument("--pdbs",dest="pdbs", nargs="+",
                    help="Provide resid mapping file generated by reduce",default=[])
parser.add_argument("--mappings",dest="mappings", nargs="+",
                    help="Provide resid mapping file generated by reduce", default=[])
parser.add_argument("--k", type=float, default=1.0,
                    help="Force constant for distance restraints")
parser.add_argument("--mode", type = str, choices = ["harmonic", "haddock", "position"], default = "harmonic",
                    help="""Specifies how the restraints in the .tbl are interpreted. The following modes are possible
                    harmonic: Defines harmonic distance restraints between two atomic selections. Format: <selection 1> <selection 2> <distance> <dminus> <dplus>
                    haddock: Same as harmonic, but enables soft-square potential and restraint cross-validation
                    position: Defines an harmonic distance restraint for each atom in a selection to a reference point in space. 
                      Format: <selection> <distance> <dminus> <dplus> <type> <x> <y> <z>                   
                      <type> can be: x, y, z, xy, xz, yz or xyz
                       The type determines which coordinate(s) of the atom are used to compute <distance>
                        e.g. "x" defines a planar region, "xy" a cylindrical region and "xyz" a spherical region                       
                      <x> <y> <z>: defines the reference point
                    """)
parser.add_argument("--softsquare", type=float, default=2.0,
                    help="HADDOCK-mode softsquare limit, indicating the maximum distance violation after which the restraint potential becomes linear")
parser.add_argument("--chance_removal", type=float, default=0.5,
                    help="chance of removal for HADDOCK-mode restraints; equal to 1/noecv, where noecv is the HADDOCK cross-validation parameter")

args = parser.parse_args()
assert len(args.pdbs) == len(args.mappings), (len(args.pdbs), len(args.mappings))
tbldata = open(args.tblfile).read()

try:
  import grako
except ImportError:
  raise ImportError("Parsing CNS files requires the Grako library")  

from tbl_grammar import tbl_grammarParser as Parser
from rmsdlib import read_pdb

pdbs = []
for p in args.pdbs:
  pdbs.append(read_pdb(p))

resmaps = []
for m in args.mappings:
  resmap = {}
  for l in open(m):
    ll = l.split()
    if not len(ll) == 2: continue
    resmap[ll[0]] = int(ll[1])
  resmaps.append(resmap)

class ATTRACTSemantics(object):

    def assign_statement_2_cross(self, ast):    
        raise NotImplementedError
        
    def assign_statement_4(self, ast):    
        raise NotImplementedError

    def assign_statement_4_cross(self, ast):    
        raise NotImplementedError

    def assign_statement_6(self, ast):    
        raise NotImplementedError

    def assign_statement_6_cross(self, ast):    
        raise NotImplementedError

    def assign_statement_pcs(self, ast):    
        raise NotImplementedError

    def assign_statement_pcs_cross(self, ast):    
        raise NotImplementedError

    def byres(self, ast):
        raise NotImplementedError

    def bygroup(self, ast):
        raise NotImplementedError

    def bondedto(self, ast):
        raise NotImplementedError

    def around(self, ast):
        raise NotImplementedError

    def saround(self, ast):
        raise NotImplementedError

    def chemical(self, ast):
        raise NotImplementedError
      
    def atom(self, ast):
        raise NotImplementedError
      
    def attribute(self, ast):
        raise NotImplementedError

    def fbox(self, ast):
        raise NotImplementedError

    def sfbox(self, ast):
        raise NotImplementedError
      
    def point(self, ast):
        raise NotImplementedError
      
    def RECALL_STORE(self, ast):
        raise NotImplementedError

    def known(self, ast):
        raise NotImplementedError
      
    def hydrogen(self, ast):
        raise NotImplementedError
      
    def all_(self, ast):
        raise NotImplementedError

    def previous(self, ast):
        raise NotImplementedError

    def tag(self, ast):
        raise NotImplementedError

    def none(self, ast):
        raise NotImplementedError
      
    def id(self, ast):
        raise NotImplementedError
      
    def _default(self, ast):
        return ast

import numpy
pdblen = [len(list(p.atoms())) for p in pdbs]
pdbcumlen = [0]
for plen in pdblen:
  pdbcumlen.append(plen + pdbcumlen[-1])  
mask0 = numpy.zeros(pdbcumlen[-1], dtype="bool")

def select_segid(segid):
  m = mask0.copy()
  segid = segid.strip()
  assert len(segid) == 1 #TODO: support for body >26
  body = ord(segid) - ord('A') + 1
  assert body <= len(pdbs), (body, len(pdbs))
  m[pdbcumlen[body-1]:pdbcumlen[body]] = 1
  #print "SEGID", segid, sum(m)
  return m
   
def select_resid(resid):
  m = mask0.copy()
  count = 0
  resid = resid.strip()
  for pnr, p in enumerate(pdbs):
    resmap = resmaps[pnr]
    resnr = None
    if resid in resmap: 
      resnr = resmap[resid]
    for a in p.atoms():
      if a.resnr == resnr: 
        m[count] = 1
      count += 1 
  #print "RESID", resid, sum(m)    
  return m

def select_resname(resname):
  m = mask0.copy()
  count = 0
  resname = resname.strip()
  for p in pdbs:
    for a in p.atoms():
      if a.resname == resname:
        m[count] = 1
      count += 1 
  return m

def select_name(name):
  m = mask0.copy()
  count = 0
  name = name.strip()
  for p in pdbs:
    for a in p.atoms():
      if a.name == name:
        m[count] = 1
      count += 1 
  return m

def evaluate_term(aa):
  term = None
  if "factor" in aa and aa["factor"] is not None:
    term = evaluate_term(aa["factor"])
  elif "term" in aa and aa["term"] is not None:
    term = evaluate_term(aa["term"])
    
  if term is not None:  
    if "or_" in aa and aa["or_"] is not None:
      aaa = aa["or_"]
      if not isinstance(aaa, list): 
        aaa = [aaa]
      for aa2 in aaa:
        term2 = evaluate_term(aa2)
        term = term | term2          
        #print "OR", sum(term)      
    if "and_" in aa and aa["and_"] is not None:      
      aaa = aa["and_"]
      if not isinstance(aaa, list): 
        aaa = [aaa]
      for aa2 in aaa:
        term2 = evaluate_term(aa2)
        term = term & term2          
        #print "AND", sum(term)      
    return term
  elif "not_" in aa and aa["not_"] is not None:
    term = evaluate_term(aa["not_"])
    return ~term
  elif "segid" in aa and aa["segid"] is not None:
    return select_segid(aa["segid"])
  elif "resid" in aa and aa["resid"] is not None:
    return select_resid(aa["resid"])
  elif "resname" in aa and aa["resname"] is not None:
    return select_resname(aa["resname"])
  elif "name" in aa and aa["name"] is not None:
    return select_name(aa["name"])
  else:
    raise Exception(aa)


#Restraint data
class Restraint(object):
  def __init__(self):
    self.maskindices = []

restraints = []
masks = []
maskhashes = []

parser = Parser(semantics=ATTRACTSemantics())
ast = parser.parse(tbldata, rule_name = "assign_statements")
assert isinstance(ast, list), type(ast)
if not len(ast):
  raise ValueError("Cannot parse TBL file")
for a in ast:
  assert isinstance(a, dict) and "assign" in a
  assign = a["assign"]
  rest = Restraint()
  rest.distance = float(a["distance"])
  rest.dminus = float(a["dminus"])
  rest.dplus = float(a["dplus"])
  
  def evaluate_assign(assi):
    mask = evaluate_term(assi)
    maskhash = hash(mask.tostring())
    #print "TERM", sum(mask)
    try:
      maskindex = maskhashes.index(maskhash)
    except ValueError:
      maskhashes.append(maskhash)
      masks.append(mask)
      maskindex = len(masks) - 1
    return maskindex
  if args.mode == "position":
    rest.xyz = a["xyz"]
    rest.typ = a["type"][0]
    maskindex = evaluate_assign(assign)
    rest.maskindices.append(maskindex)
  else: #"harmonic", "haddock"
    assert len(assign) == 2  
    for aa in assign:
      assert isinstance(aa, dict), "Malformed assign statement"
      maskindex = evaluate_assign(aa)    
      rest.maskindices.append(maskindex)
  restraints.append(rest)
  
for mnr, m in enumerate(masks):
  print "selection%d" % (mnr+1),  
  print sum(m),
  for mm in numpy.where(m>0)[0]: print mm+1,
  print
print
for rest in restraints:
  mindist = rest.distance - rest.dminus
  maxdist = rest.distance + rest.dplus
  if args.mode == "position":
    x, y, z = rest.xyz.x, rest.xyz.y, rest.xyz.z
    r = "selection%d 7 %s %s %s %s %s %s %s" % (rest.maskindices[0]+1, mindist, maxdist, args.k, rest.typ, x, y, z)
    print r
  else: # "harmonic", "haddock"
    r = "selection%d selection%d" % (rest.maskindices[0]+1, rest.maskindices[1]+1)  
    if args.mode == "haddock": 
      if mindist > 0:
        raise Exception("HADDOCK restraints with minimum distance > 0 are not supported")
      print r + " 2 %s %s %s %s" % (maxdist, args.k, args.softsquare, args.chance_removal)
    else: #harmonic
      if mindist == maxdist:
        print r + " 4 %s %s" % (maxdist, args.k)
      else:
        print r + " 1 %s %s" % (maxdist, args.k)
        if mindist > 0:
          print r + " 3 %s %s" % (mindist, args.k)
      
  
