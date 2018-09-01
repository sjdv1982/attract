class Residue(object):
    topname = None #name of topology where this residue has been defined, e.g. "oplsx" or "dna-rna"
    def __init__(self, name, atomorder, atoms, bonds, topname=None):
        self.topname = topname
        self.name = name
        self.atomorder = atomorder
        self.atoms = atoms
        self.bonds = bonds
        self.patches = []

    def add_atom(self, atom):
        if atom.name not in self.atomorder:
            self.atomorder.append(atom.name)
        self.atoms[atom.name] = atom

    def add_bond(self, bond):
        atom1, atom2 = bond
        self.bonds.append((atom1, atom2))

    def delete_atom(self, name):
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

    def modify_atom(self, atom):
        assert atom["name"] in self.atoms, atom["name"]
        self.atoms[atom["name"]].update(atom)

    def patch(self, patch):
        name, commands = patch
        if name in self.patches:
            return
        self.patches.append(name)
        for mode, obj in commands:
            if mode == "add_atom":
                self.add_atom(obj)
            elif mode == "add_bond":
                self.add_bond(obj)
            elif mode == "modify_atom":
                self.modify_atom(obj)
            elif mode == "delete_atom":
                self.delete_atom(obj)
            else:
                raise ValueError(mode)

def load(d):
    residues = {}
    for resname, resdata in d["residues"].items():
        res = Residue(resname, **resdata)
        residues[resname] = res
    patches = d["patches"]
    return residues, patches

def merge(forcefields):
    residues, patches = {}, {}
    for new_residues, new_patches in forcefields:
        for resname in new_residues:
            new_res = new_residues[resname]
            if resname in residues:
                res = residues[resname]
                topname = res.topname
                new_topname = new_res.topname
                if topname == new_topname:
                    msg = "both have topology name '%s'" % topname
                else:
                    msg = "topology names: '%s' and '%s'" % (topname, new_topname)
                raise ValueError("Duplicate residue %s, %s" % (resname, msg))
            residues[resname] = new_res
        for patchname in new_patches:
            new_patch = new_patches[patchname]
            if patchname in patches:
                patch = patches[patchname]
                topname = patch["topname"]
                new_topname = new_patch["topname"]
                if topname == new_topname:
                    msg = "both have topology name '%s'" % topname
                else:
                    msg = "topology names: '%s' and '%s'" % (topname, new_topname)
                raise ValueError("Duplicate patch %s, %s" % (patchname, msg))
            patches[patchname] = new_patch
    patches = {k:(k, v) for k,v in patches.items()}
    return residues, patches
