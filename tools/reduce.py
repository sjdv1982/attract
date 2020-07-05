"""
Reduce converts a PDB into a modified PDB format, where the atom type for each atom is indicated
This reduced PDB is what is understood by ATTRACT.
These atoms are actually pseudo-atoms (beads), representing the average position of one or more real atoms.
The pseudo-atoms definitions are in reduce.param, see that file for more details.

DNA bases are DA, DC, DG, DT; and RNA bases are RA, RC, RG, RU. One-letter and three-letter nucleic acid\
codes can be  automatically interpreted as DNA or RNA with the --dna and --rna options.

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

def read_filelist(filelist):
    ret = []
    for l in open(filelist):
        l = l.strip()
        if not len(l): continue
        assert len(l.split()) == 1, (filelist, l)
        ret.append(l)
    return ret

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


def print_res(mapping=True):
    if not len(state.rescoor): return
    state.rescounter += 1
    if mapping and state.mapf is not None:
        print(state.res[1:].strip(), state.rescounter, file=state.mapf)
    for l in state.beadgroups[state.resname]:
        if (l[0], l[1]) not in state.rescoor:
            continue
        c = state.rescoor[(l[0], l[1])]
        x, y, z = c[1]/c[0], c[2]/c[0], c[3]/c[0]
        state.atomcounter += 1
        atomname = l[1]
        if len(atomname) < 4:
                atomname = " " + atomname + "   "[len(atomname):]
        line = (state.atomcounter, atomname, state.resname, state.chain, state.rescounter, x, y, z, l[0], l[3], 1.0)
        print("ATOM%7d %4s %3s %s%4d    %8.3f%8.3f%8.3f%5d%8.3f 0%5.2f" % line, file=state.outp)
    state.rescoor = {}

def run(pdbdata,mapping=True):
    state.res = None
    state.resname = None
    state.rescounter = state.startres-1
    state.atomcounter = state.startatom-1
    state.rescoor = {}


    for l in pdbdata.splitlines():
        if not l.startswith("ATOM"): continue
        cres = l[21:26]
        if cres != state.res:
            print_res(mapping)
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
        for bead in bgres:
            for at in bead[2]:
                if atom != at: continue
                beadname = bead[0], bead[1]
                if beadname not in state.rescoor:
                    state.rescoor[beadname] = [0, 0.0, 0.0, 0.0]
                c = state.rescoor[beadname]
                c[0] += 1
                c[1] += x
                c[2] += y
                c[3] += z
                break

if __name__ == "__main__":
    import argparse
    global args

    parser =argparse.ArgumentParser(description=__doc__,
                                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdb",help="PDB file to reduce")
    parser.add_argument("output",help="reduced output PDB file", nargs="?")

    parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
    parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")
    parser.add_argument("--chain", help="Set the chain in the output PDB", default=" ")
    parser.add_argument("--startres", help="Set residue number of the first residue, default 1", type=int,default=1)
    parser.add_argument("--startatom", help="Set atom number of the first atom, default 1", type=int,default=1)
    parser.add_argument("--batch", help="run reduce in batch mode. Input and output must be two lists of PDBs", action="store_true")

    args = parser.parse_args()
    state.startres = args.startres
    state.startatom = args.startatom
    state.mapf = sys.stdout

    if args.dna and args.rna:
        raise ValueError("Options --dna and --rna are mutually exclusive")

    if args.batch and args.output is None:
        raise ValueError("--batch requires a file list as output argument")

    state.dna = args.dna
    state.rna = args.rna

    assert len(args.chain) == 1, args.chain
    state.chain = args.chain

    reduce_param = os.path.split(os.path.abspath(__file__))[0] + os.sep + "reduce-beadgroups.txt"
    beadgroups_data = open(reduce_param).read()
    state.beadgroups = read_beadgroups(beadgroups_data)

    if args.batch:
        infiles = read_filelist(args.pdb)
        for f in infiles:
            assert os.path.exists(f), f
        outfiles = read_filelist(args.output)
        for pdb, outfile in zip(infiles, outfiles):
            state.outp = open(outfile, "w")
            pdbdata = open(pdb).read()
            run(pdbdata,mapping=False)
            print_res(mapping=False)
            state.outp.close()
    else:
        pdb = args.pdb
        assert os.path.exists(pdb)
        if args.output is None:
            args.output = os.path.splitext(pdb)[0] + "r.pdb"
        state.outp = open(args.output, "w")
        pdbdata = open(pdb).read()
        run(pdbdata)
        print_res()

if __name__ == "transformer":
    """For use from Seamless
    Input pins:
    - pdbdata
    - molecule: "dna", "rna", "protein"
    - chain (normally " ")
    - starting_residue (normally 1)
    - starting_atomindex (normally 1)
    - beadgroups
    """
    from io import StringIO
    state = State()
    state.dna = (molecule == "dna")
    state.rna = (molecule == "rna")
    state.chain = chain
    state.startres = starting_residue
    state.startatom = starting_atomindex
    state.beadgroups = read_beadgroups(beadgroups)


    state.outp = StringIO()
    state.mapf = StringIO()

    run(pdbdata)
    print_res()

    result = {}
    result["pdbdata"] = state.outp.getvalue()
    mapping0 = state.mapf.getvalue()
    mapping = {}
    for line in mapping0.splitlines():
        res_from, res_to = line.split()
        res_to = int(res_to)
        mapping[res_from] = res_to
    result["mapping"] = mapping
