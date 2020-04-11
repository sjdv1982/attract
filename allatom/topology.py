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
        name = str(atom["name"])
        if name not in self.atomorder:
            self.atomorder.append(name)
        self.atoms[name] = atom

    def add_bond(self, bond):
        atom1, atom2 = bond
        self.bonds.append((str(atom1), str(atom2)))

    def delete_atom(self, name):
        name = str(name)
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
        name = str(atom["name"])
        assert name in self.atoms, name
        self.atoms[name].update(atom)

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

    def get_angles(self):
        connectivity = {}
        angles = []
        for atom1, atom2 in self.bonds:
            try:
              a1 = self.atomorder.index(atom1)
              a2 = self.atomorder.index(atom2)
            except ValueError:
              pass
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
