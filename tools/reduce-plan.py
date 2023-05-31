"""
Version of reduce.py that outputs not a reduced PDB but a "plan": for each bead, a list of atom indices is printed.
This list indicates which atoms of the input PDB are to be averaged to compute the reduced coordinate
This list is space-separated, one list per line, indices start from 1.
This plan is used by reduce-npy.pdb to reduce all-atom Numpy arrays of coordinates

Copyright Sjoerd de Vries, Creative Commons CC-BY-SA license
"""

from __future__ import print_function

import sys, os

class State():
    pass
state = State()

#Mapping of nucleic-acid codes to DNA/RNA
mapnuc = {
    "A": ["DA", "RA"],
    "A3": ["DA", "RA"],
    "A5": ["DA", "RA"],
    "DA3": ["DA", None],
    "DA5": ["DA", None],
    "ADE": ["DA", "RA"],
    "C": ["DC", "RC"],
    "C3": ["DC", "RC"],
    "C5": ["DC", "RC"],
    "DC3": ["DC", None],
    "DC5": ["DC", None],
    "CYT": ["DC", "RC"],
    "G": ["DG", "RG"],
    "G3": ["DG", "RG"],
    "G5": ["DG", "RG"],
    "DG3": ["DG", None],
    "DG5": ["DG", None],
    "GUA": ["DG", "RG"],
    "T": ["DT", None],
    "T3": ["DT", None],
    "T5": ["DT", None],
    "DT3": ["DT", None],
    "DT5": ["DT", None],
    "THY": ["DT", None],
    "U": [None, "RU"],
    "U3": [None, "RU"],
    "U5": [None, "RU"],
    "URA": [None, "RU"],
    "URI": [None, "RU"],
}

def read_beadgroups(beadgroups_data):
    bg = {}
    aa = None
    for l in beadgroups_data.splitlines():
        pound = l.find("#")
        if pound > -1: l = l[:pound]
        l = l.strip()
        if not len(l): continue
        ll = l.split()
        if len(ll) == 1:
            aa = ll[0]
            assert len(aa) <= 3, l
            bg[aa] = []
        else:
            assert aa is not None
            try:
                atomtype = int(ll[0])
            except ValueError:
                raise ValueError(l)
            atoms = ll[2:]
            charge = 0.0
            try:
                charge = float(atoms[-1])
                atoms = atoms[:-1]
            except ValueError:
                pass
            bg[aa].append( (int(ll[0]), ll[1], set(atoms), charge) )
    return bg


def print_res():
    if not len(state.rescoor): return
    for l in state.beadgroups[state.resname]:
        if (l[0], l[1]) not in state.rescoor:
            raise ValueError("Missing {}, {} (bead code{})".format(state.resname, l[1], l[0]))
        indices = state.rescoor[(l[0], l[1])]
        for ind in indices:
            print(ind, end=" ", file=state.outp)
        print(file=state.outp)
    state.rescoor = {}

def run(pdbdata):
    state.res = None
    state.resname = None
    state.rescoor = {}


    atom_index = 0
    for l in pdbdata.splitlines():
        if not l.startswith("ATOM"): continue
        cres = l[21:26]
        if cres != state.res:
            print_res()
            state.res = cres
            state.resname = l[17:20].strip()
            if state.resname in mapnuc:
                if state.dna:
                    state.resname = mapnuc[state.resname][0]
                elif state.rna:
                    state.resname = mapnuc[state.resname][1]
                else:
                    msg = "PDB contains a nucleic acid named \"%s\", but it could be either RNA or DNA. Please specify the --dna or --rna option"
                    raise ValueError(msg % state.resname)
                if state.resname is None:
                    if state.dna: na = "DNA"
                    if state.rna: na = "RNA"
                    raise ValueError("'%s' can't be %s" % (l[17:20].strip(), na))

            assert state.resname in state.beadgroups, l
            bgres = state.beadgroups[state.resname]
        try:
            atom = l[12:16].strip()
            x = float(l[30:38])
            y = float(l[38:46])
            z = float(l[46:54])
        except ValueError:
            continue
        atom_index += 1
        for bead in bgres:
            for at in bead[2]:
                if atom != at: continue
                beadname = bead[0], bead[1]
                if beadname not in state.rescoor:
                    state.rescoor[beadname] = []
                state.rescoor[beadname].append(atom_index)
                break

if __name__ == "__main__":
    import argparse
    global args

    parser =argparse.ArgumentParser(description=__doc__,
                                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdb",help="PDB file to reduce")
    parser.add_argument("output",help="output plan file", nargs="?")

    parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
    parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")


    args = parser.parse_args()
    state.mapf = sys.stdout

    if args.dna and args.rna:
        raise ValueError("Options --dna and --rna are mutually exclusive")

    state.dna = args.dna
    state.rna = args.rna


    reduce_param = os.path.split(os.path.abspath(__file__))[0] + os.sep + "reduce-beadgroups.txt"
    beadgroups_data = open(reduce_param).read()
    state.beadgroups = read_beadgroups(beadgroups_data)

    pdb = args.pdb
    assert os.path.exists(pdb)
    if args.output is None:
        state.outp = sys.stdout
    else:
        state.outp = open(args.output, "w")
    pdbdata = open(pdb).read()
    run(pdbdata)
    print_res()

if __name__ == "transformer":
    """For use from Seamless
    Input pins:
    - pdbdata
    - molecule: "dna", "rna", "protein"
    - beadgroups
    """
    from io import StringIO
    state = State()
    state.dna = (molecule == "dna")
    state.rna = (molecule == "rna")
    state.beadgroups = read_beadgroups(beadgroups)


    state.outp = StringIO()

    run(pdbdata)
    print_res()

    result = state.outp.getvalue()
