#CNS parameter file parser, May 2011, Sjoerd de Vries
#Reads in residue and residue patch (presidue) statements
#Currently, all but atom-related statements are ignored
#Note: +C and -C modifiers are ignored!

import json

allspecies = ("atom", "angle", "bond", "improper","impr","dono", "donor", "acce", "acceptor", "dihed", "dihedral", "group", "grou")

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
    def serialize(self):
        result = {}
        if self.name is not None:
            result["name"] = self.name
        if self.type is not None:
            result["type"] = self.type
        if self.charge is not None:
            result["charge"] = self.charge
        return result

class Residue(object):
    topname = None #name of topology where this residue has been defined, e.g. "oplsx" or "dna-rna"
    def __init__(self, name):
        self.name = name
        self.atomorder = []
        self.atoms = {}
        self.bonds = []

    def _identify(self, species):
        assert species in allspecies, species
        if species == "atom":
            order = self.atomorder
            data = self.atoms
        elif species == "bond":
            order = [] #dummy
            data = self.bonds
        else:
            return None
        return data, order

    def add(self, species,tokens):
        id = self._identify(species)
        if id is None: return
        data, order = id
        if species == "atom":
            assert tokens[-1] == "end", tokens[-1]
            a = Atom.parse(tokens[:-1])
            order.append(a.name)
            data[a.name] = a
        else:
            a = tokens
            if len(a) != 2 and species == "bond":
                assert len(a) >= 5 and a[2] == "bond", a
                data.append(a[:2])
                self.add("bond", a[3:])
            else:
                if species == "bond":
                    a = [aa.lstrip("-").lstrip("+") for aa in a]
                data.append(a)

    def delete(self, species, name):
        assert species == "atom", species
        try:
            del self.atoms[name]
            self.atomorder.remove(name)
        except KeyError:
            pass
        bonds2 = []
        for b in self.bonds:
            if b[0] != name and b[1] != name: bonds2.append(b)
        #TODO: delete dihedrals (once we read them)
        self.bonds = bonds2

    def modify(self, species, tokens):
        assert species == "atom", species
        assert tokens[-1] == "end", tokens
        if len(tokens) == 2:
            return
        a = Atom.partparse(tokens[:-1])
        assert a.name in self.atoms, a.name
        self.atoms[a.name].modify(a)

    def copy(self):
        ret = type(self)(self.name)
        ret.topname = self.topname
        ret.atomorder = list(self.atomorder)
        for aname in self.atoms:
            ret.atoms[aname] = self.atoms[aname].copy()
        ret.bonds[:] = self.bonds
        return ret

    def serialize(self):
        atoms = {k:v.serialize() for k,v in self.atoms.items()}
        return {
            "topname": self.topname,
            #"name": self.name, #in parent
            "atomorder": self.atomorder,
            "atoms": atoms,
            "bonds": self.bonds,
        }

    def get_angles(self):
        """Not used by ATTRACT..."""
        connectivity = {}
        angles = []
        for atom1, atom2 in self.bonds:
            a1 = self.atomorder.index(atom1)
            a2 = self.atomorder.index(atom2)
            if atom1 not in connectivity: connectivity[atom1] = set()
            if atom2 not in connectivity: connectivity[atom2] = set()
            con1 = connectivity[atom1]
            con2 = connectivity[atom2]
            con1.add(a2)
            con2.add(a1)
        for atom in connectivity.keys():
            con = connectivity[atom]
            for a1 in con:
                atom1 = self.atomorder[a1]
                for a2 in con:
                    if a2 <= a1: continue
                    atom2 = self.atomorder[a2]
                    angles.append((atom1, atom, atom2))
        return angles



class Patch(object):
    topname = None #name of topology where this residue has been defined, e.g. "oplsx" or "dna-rna"
    def __init__(self, name):
        self.name = name
        self.commands = []
    def add_command(self,mode,tokens):
        species = tokens[0]
        assert species in allspecies
        if species not in ("bond", "atom"): return
        if mode == "dele": mode = "delete"
        if mode == "modi": mode = "modify"
        if mode == "add":
            if species == "atom":
                atom = Atom.parse(tokens[1:])
                self.commands.append(("add_atom", atom.serialize()))
            else:
                assert len(tokens) == 3, tokens
                bond = [aa.lstrip("-").lstrip("+") for aa in tokens[1:]]
                self.commands.append(("add_bond", bond))
        elif mode == "modify":
            assert tokens[-1] == "end"
            if len(tokens) == 3:
                return
            name = tokens[1].lstrip("-").lstrip("+")
            #kludge to deal with 2-residue patches (e.g. DISU)
            #modifications for group 1 are taken, for group 2 are ignored
            if name.startswith("1"): name = name[1:]
            if name.startswith("2"): return
            if species != "atom":
                return
            atom = Atom.partparse([name] + tokens[2:-1])
            self.commands.append(("modify_atom", atom.serialize()))
        elif mode == "delete":
            assert tokens[-1] == "end", tokens
            assert len(tokens) == 3
            name = tokens[1].lstrip("-").lstrip("+")
            #kludge to deal with 2-residue patches (e.g. DISU)
            #modifications for group 1 are taken, for group 2 are ignored
            if name.startswith("1"): name = name[1:]
            if name.startswith("2"): return
            if species != "atom": return
            self.commands.append(("delete_atom", name))
        else:
            raise ValueError("unknown mode: %s", mode)
    def serialize(self):
        return self.commands
residues = {}
patches = {}
def parse_stream(stream,topname=None):
    mode = 0
    linenr = 0
    for line in stream:
        line = line.replace("{","").replace("}","")
        linenr += 1
        l = line.rstrip("\n").rstrip("\r")
        tokens = tokenize(l)
        if len(tokens) == 0: continue
        if mode == 0:
            if tokens[0] in ("residue", "resi", "resid"):
                name = tokens[1]
                res = Residue(name)
                res.topname = topname
                residues[name] = res
                mode = 1; continue
            elif tokens[0] in ("presidue", "pres", "presi","presid"):
                name = tokens[1]
                res = Patch(name)
                res.topname = topname
                patches[name] = res
                mode = 2; continue
            else: continue
        elif mode == 1:
            if tokens[0] != "group": continue
            mode = 2; continue
        if tokens[0] == "end":
            mode = 0; continue
        #print linenr, tokens
        typename = res.__class__.__name__
        if typename == "Residue":
            res.add(tokens[0],tokens[1:])
        elif typename == "Patch":
            if tokens[0] == "group": continue
            try:
                res.add_command(tokens[0], tokens[1:])
            except:
                print "Line %d: %s" % (linenr, line)
                raise

        #print res.__class__.__name__, name, tokens

if __name__ == "__main__":
    import sys, os
    topname = sys.argv[1]
    topfiles = sys.argv[2:]
    for topfile in topfiles:
        stream = open(topfile).readlines()
        parse_stream(stream, topname=topname)
    residues = {k.upper(): v.serialize() for k,v in residues.items()}
    patches = {k.upper(): v.serialize() for k,v in patches.items()}
    result = {
        "residues": residues,
        "patches": patches,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
