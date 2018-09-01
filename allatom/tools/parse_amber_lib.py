import json
import parmed
import sys

residues = {}
topname = sys.argv[1]
libfiles = sys.argv[2:]
for libfile in libfiles:
    lib = parmed.amber.offlib.AmberOFFLibrary.parse(libfile)
    for resname in lib:
        res = lib[resname]
        atoms = {}
        atomorder = []
        r = {
            "topname": topname,
            "atoms": atoms,
            "atomorder": atomorder,
            "bonds": [], #not implemented yet
        }
        for atom in res.atoms:
            name = atom.name
            a = {
                "name": name,
                "type": atom.type.upper(),
                "charge": atom.charge,
            }
            atoms[name] = a
            atomorder.append(name)
        residues[resname] = r

result = {
    "residues": residues,
    "patches": [],
}
print(json.dumps(result, indent=2, sort_keys=True))
