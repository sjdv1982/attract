#CNS parameter file parser, May 2011, Sjoerd de Vries
#Reads in residue and residue patch (presidue) statements
#Currently, all but atom-related statements are ignored
#Note: +C and -C modifiers are ignored! 

allspecies = ("atom", "angle", "bond", "improper","impr","dono", "donor", "acce", "acceptor", "dihed", "dihedral")

def tokenize(s):
  s = s.lower()
  ss = s.split()
  r = []
  for t in ss:
    pos = t.find("=")
    while pos > -1:
      r.append(t[:pos])
      r.append("=")
      t = t[pos+1:]
      pos = t.find("=")
    r.append(t)
  ret = []
  for rr in r:
    if rr.startswith("!"): break
    if len(rr) == 0: continue
    ret.append(rr)
  return ret

class Atom(object):
  def __init__(self, name, type, charge):
    self.name = name.lstrip("-").lstrip("+")
    self.type = type
    self.charge = charge
  def copy(self):
    return type(self)(self.name, self.type, self.charge)
  def modify(self, a):
    if a.name is not None: self.name = a.name
    if a.type is not None: self.type = a.type
    if a.charge is not None: self.charge = a.charge    
  @staticmethod
  def _parse(tokens):
    name = tokens[0]
    d = {"type":None,"charge":None}
    assert tokens[1] in ("type","charge"), tokens
    assert tokens[2] == "=", tokens
    d[tokens[1]] = tokens[3]
    if len(tokens) > 4:
      assert tokens[4] in ("type","charge"), tokens
      assert tokens[5] == "=", tokens
      d[tokens[4]] = tokens[6]
    if d["charge"] is not None: d["charge"] = float(d["charge"])
    return name, d["type"], d["charge"]
  @classmethod
  def parse(cls, tokens):
    name, type, charge = cls._parse(tokens)
    assert type is not None and charge is not None, (type, charge)
    return cls(name, type, charge)
  @classmethod
  def partparse(cls, tokens):
    name, type, charge = cls._parse(tokens)
    return cls(name, type, charge)
  def __str__(self):
    return "atom %s type=%s charge=%.6f" % (self.name, self.type, self.charge)
        
class Residue(object):
  def __init__(self, name):
    self.name = name
    self.atoms = {}
  def add(self, species,tokens):
    assert species in allspecies, species
    if species != "atom": return
    assert tokens[-1] == "end", tokens[-1]
    a = Atom.parse(tokens[:-1])
    self.atoms[a.name] = a
  def delete(self, species, name):
    assert species in allspecies
    if species != "atom": return
    assert name in self.atoms, name
    del self.atoms[name]
  def modify(self, species, a):
    assert species in allspecies
    if species != "atom": return
    assert a.name in self.atoms, a.name
    self.atoms[a.name].modify(a)
  def copy(self):
    ret = type(self)(self.name)
    for aname in self.atoms:
      ret.atoms[aname] = self.atoms[aname].copy()
    return ret
  def patch(self, pres):
    for mode, species, obj in pres.commands:
      if species != "atom": continue
      if mode == "add":
        self.atoms[obj.name] = obj
      elif mode == "modify":
        self.modify(species, obj)
      elif mode == "delete":
        self.delete(species, obj)
      else: raise ValueError(mode)
    

class PResidue(object):
  def __init__(self, name):
    self.name = name
    self.commands = []
  def add_command(self,mode,tokens):
    species = tokens[0]
    assert species in allspecies
    if species != "atom": return
    if mode == "dele": mode = "delete"
    if mode == "modi": mode = "modify"
    if mode == "add":
      a = Atom.parse(tokens[1:])
      self.commands.append((mode, species, a))
    elif mode == "modify":
      assert tokens[-1] == "end"
      a = Atom.partparse(tokens[1:-1])
      self.commands.append((mode, species, a))
    elif mode == "delete":
      assert tokens[-1] == "end", tokens
      assert len(tokens) == 3
      name = tokens[1].lstrip("-").lstrip("+")
      self.commands.append((mode, species, name))
    else:
      raise ValueError("unknown mode: %s", mode)

def parse_stream(stream): 
  residues = {}
  presidues = {}
  mode = 0
  linenr = 0
  for line in stream:
    linenr += 1
    l = line.rstrip("\n").rstrip("\r")
    tokens = tokenize(l)
    if len(tokens) == 0: continue
    if mode == 0:
      if tokens[0] in ("residue", "resi", "resid"): 
        name = tokens[1]
        res = Residue(name)
        residues[name] = res
	mode = 1; continue
      elif tokens[0] in ("presidue", "pres", "presi","presid"):
        name = tokens[1]
        res = PResidue(name)
	presidues[name] = res
        mode = 2; continue      
      else: continue
    elif mode == 1:
      if tokens[0] != "group": continue
      mode = 2; continue
    if tokens[0] == "end":
      mode = 0; continue
    print linenr, tokens
    typename = res.__class__.__name__
    if typename == "Residue":
      res.add(tokens[0],tokens[1:])
    elif typename == "PResidue":
      if tokens[0] == "group": continue
      res.add_command(tokens[0], tokens[1:])
      
    #print res.__class__.__name__, name, tokens
  return residues, presidues
    
    
if __name__ == "__main__":
  residues, presidues = parse_stream(open("topallhdg5.3.pro"))      
  print residues.keys()    
  a = residues["ala"].copy()
  for atom in a.atoms: print a.atoms[atom]
  print  
  a.patch(presidues["cter"])
  for atom in a.atoms: print a.atoms[atom]
