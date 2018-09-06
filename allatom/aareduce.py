#!/usr/bin/env python2
"""
Converts a PDB into a modified PDB format, where the atom type for each atom is indicated
This reduced PDB is what is understood by ATTRACT.
These atoms are parsed from the trans files and topology files derived from the OPLS forcefield

For DNA and RNA nucleotides, "dna-rna-top.json" is read in
These nucleotides, by default, contain a phosphate before the sugar
Since this phosphate is often missing for the first nucleotide, the patch "5ter" is applied by default
To prevent this, use "--patch 1 None"; this will keep the phosphate, but not the O5T
To keep the O5T as well, use the --termini or --nter option.
"""
from __future__ import print_function
import sys, os, json
import pdbcomplete
import topology
from pdbcomplete import run_pdb2pqr, pdbfix, update_patches, pdb_lastresort
from copy import deepcopy

has_argparse = False
try:
    import argparse
    has_argparse = True
except ImportError:
    import optparse  #Python 2.6

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
mapnucrev = {
  "DA":"A",
  "RA":"A",
  "DC":"C",
  "RC":"C",
  "DG":"G",
  "RG":"G",
  "DT":"T",
  "RU":"U",
}
class PDBres:
    def __init__(self, chain, resid, resname, topology):
        self.chain = chain
        self.resid = resid
        self.resname = resname
        self.coords = {}
        self.chainfirst = False
        self.chainlast = False
        self.nter = False
        self.cter = False
        self.topology = topology

code_to_type = {}
def parse_transfile(transfile, topname):
    for l in open(transfile):
        ll = l.split()
        type = int(ll[0])
        for code in ll[3:]:
            if code.startswith("#"): break
            assert (code, topname) not in code_to_type, code
            code_to_type[code, topname] = type
            code_to_type[str(code), str(topname)] = type

mutations = {}

def read_filelist(filelist):
    ret = []
    for l in open(filelist):
        l = l.strip()
        if not len(l): continue
        assert len(l.split()) == 1, (filelist, l)
        ret.append(l)
    return ret

def read_pdb(pdblines, pdbname, add_termini=False,modbase=False,modres=False):
    repl = (
      ("H","HN"),
      ("HT1","HN"),
      ("OP1","O1P"),
      ("OP2","O2P"),
      ("H1","HN"),
      ("OP3","O5T"),
      ("HO5'","H5T"),
      ("HO3'","H3T"),
    )
    pdbres = []
    curr_res = None
    atomlines = []

    if (modbase or modres):
        res0 = {}
        pdblines = list(pdblines)
        for l in pdblines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                resid = l[22:27]
                if resid not in res0: res0[resid] = set()
                atomcode = l[12:16].strip()
                res0[resid].add(atomcode)
        res_ok = set()
        for r in res0:
            ratoms = res0[r]
            if modres and "CA" in ratoms and "C" in ratoms and "N" in ratoms:
                res_ok.add(r)
            elif modbase:
                cp = 0
                for n in range(5):
                    if ("C%d'" % n) in ratoms: cp += 1
                if cp >= 3:
                    res_ok.add(r)
        for l in pdblines:
            if l.startswith("ATOM"):
                atomlines.append(l)
            elif l.startswith("HETATM"):
                resid = l[22:27]
                if resid in res_ok:
                    atomlines.append(l)
    else:
        atomlines = [l for l in pdblines if l.startswith("ATOM")]

    if len(atomlines) == 0:
        raise ValueError("PDB '%s' contains no ATOM records" % pdbname )

    for l in atomlines:
        atomcode = l[12:16].strip()
        if l[16] not in (" ", "A"): continue #only keep the first of alternative conformations
        if l[30:38] == " XXXXXXX": continue #missing atom from --manual mode
        resname = l[17:20].strip()
        if resname in ["HIE", "HIP", "HSD", "HSE"]: resname="HIS"
        if resname=="HYP" and atomcode=='OD1':atomcode='OG1'
        if resname=="HYP" and atomcode=='HD1':atomcode='HG1'
        if resname in mutations: resname = mutations[resname]
        if resname in mapnuc:
            if args.dna:
                resname = mapnuc[resname][0]
            elif args.rna:
                resname = mapnuc[resname][1]
            else:
                raise ValueError("PDB contains a nucleic acid named \"%s\", but it could be either RNA or DNA. Please specify the --dna or --rna option" % resname)
            if resname is None:
                if args.dna: na = "DNA"
                if args.rna: na = "RNA"
                raise ValueError("'%s' can't be %s" % (l[17:20].strip(), na))
        chain = l[21]
        resid = l[22:27]
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        newres = False
        nter = False
        chainfirst = False
        if curr_res is None:
            newres = True
            chainfirst = True
            if add_termini: nter = True
        elif chain != curr_res.chain:
            newres = True
            chainfirst = True
            curr_res.chainlast = True
            if add_termini:
                nter = True
                curr_res.cter = True
        elif resid != curr_res.resid or resname != curr_res.resname:
            newres = True
        if newres:
            try:
                if resname is None: raise KeyError
                topr = deepcopy(top_residues[resname])
            except KeyError:
                raise KeyError("Residue type %s not known by the topology file" % resname)
            curr_res = PDBres(chain, resid, resname, topr)
            if chainfirst: curr_res.chainfirst = True
            if nter: curr_res.nter = True
            pdbres.append(curr_res)
        curr_res.coords[atomcode] = (x,y,z)
        for pin, pout in repl:
            if atomcode != pin: continue
            curr_res.coords[pout] = (x,y,z)
    if curr_res is not None:
        curr_res.chainlast = True
        if add_termini:
            curr_res.cter = True
    return pdbres

def termini_pdb(pdbres, nter, cter):
    xter = nter, cter
    for n in range(2):
        ter = xter[n]
        for resnr in ter:
            r = [res for res in pdbres if res.resid == resnr]
            if len(r) == 0:
                raise ValueError("Cannot find residue %d" % resnr)
            elif len(r) > 1:
                raise ValueError("Multiple residues %d" % resnr)
            res = r[0]
            if n == 0: res.nter = True
            else: res.cter = True

def patch_pdb(pdbres, patches):
    for res in pdbres:
        if res.resid in patches:
            for p in patches[res.resid]:
                if p is None: continue
                res.topology.patch(top_patches[p.upper()])
        elif len(pdbres) > 1 and "ca" in res.topology.atomorder: #protein
            if res.nter:
                if res.resname == "PRO":
                    res.topology.patch(top_patches["PROP"])
                else:
                    res.topology.patch(top_patches["NTER"])
            if res.cter:
                res.topology.patch(top_patches["CTER2"])
        elif len(pdbres) > 1 and "p" in res.topology.atomorder: #DNA/RNA
            if res.chainfirst:
                if res.nter:
                    res.topology.patch(top_patches["5PHO"])
                else:
                    res.topology.patch(top_patches["5TER"])
            if res.chainlast:
                res.topology.patch(top_patches["3TER"])

def check_pdb(pdbres, heavy=False):
    for res in pdbres:
        top = res.topology
        for a in top.atomorder:
            atom = top.atoms[a]
            if a.lower().startswith("h"):
                if heavy: continue
                if atom["charge"] == 0: continue
            aa = a.upper()
            if aa.strip() not in res.coords:
                raise ValueError('Missing coordinates for atom "%s" in residue %s %s%s' % (aa.strip(), res.resname, res.chain, res.resid))

def write_pdb(pdbres, chain, heavy = False, one_letter_na = False):
    pdblines = []
    mapping = []
    atomcounter = args.startatom
    rescounter = args.startres
    for res in pdbres:
        top = res.topology
        for a in top.atomorder:
            atom = top.atoms[a]
            if a.lower().startswith("h"):
                if heavy: continue
                if atom["charge"] == 0: continue
            aa = a.upper()
            x = " XXXXXXX"
            y = x; z = x
            if aa.strip() in res.coords:
                x,y,z = ("%8.3f" % v for v in res.coords[aa.strip()])
            xyz = x + y + z
            type = code_to_type[atom["type"].upper(), top.topname]
            a0 = aa
            if len(a0) < 4:
                a0 = " " + a0 + "   "[len(a0):]
            resname = res.resname
            atomname = a0
            if one_letter_na and resname in mapnucrev:
                resname = mapnucrev[resname]
            pdblines.append("ATOM%7d %4s %3s %s%4d    %s %4d %7.3f 0 1.00" % \
              (atomcounter, atomname, resname, chain, rescounter, xyz, type, atom["charge"]))
            atomcounter += 1
        mapping.append((res.resid, rescounter))
        rescounter += 1
    return pdblines, mapping

def set_reference(pdbres, pdbreferes):
    if len(pdbres) != len(pdbreferes):
        raise ValueError("PDB and reference do not have the same number of residues, %d vs %s" % (len(pdbres), len(pdbreferes)))
    for n in range(len(pdbres)):
        pdbr, refr = pdbres[n], pdbreferes[n]
        if pdbr.resname != refr.resname:
            rsid = pdbr.resid
            if refr.resid != pdbr.resid: rsid = "%s(%s)" % (pdbr.resid, refr.resid)
            raise ValueError("PDB and reference are different at resid %s: %s vs %s" % (rsid, pdbr.resname, refr.resname))
        pdbr.nter = refr.nter
        pdbr.cter = refr.cter
        pdbr.topology = refr.topology

currdir = os.path.abspath(os.path.split(__file__)[0])
topfiles = [ currdir + "/oplsx-top.json",
              currdir + "/dna-rna-top.json",
            ]
transfiles = [currdir + "/oplsx.trans",
              currdir + "/dna-rna.trans"
             ]

if has_argparse:
    parser = argparse.ArgumentParser(description=__doc__,
                              formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pdb",help="PDB file to reduce")
    parser.add_argument("output",help="all-atom reduced output PDB file", nargs="?")
else:
    parser = optparse.OptionParser()
    parser.add_argument = parser.add_option

parser.add_argument("--heavy",help="Ignore all hydrogens", action="store_true")
parser.add_argument("--refe", "--reference",help="Analyze the hydrogens of a reference file to determine histidine/cysteine states")
parser.add_argument("--autorefe",help="Analyze the hydrogens of the input PDB to determine histidine/cysteine states", action="store_true")
parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")
parser.add_argument("--nalib",help="Use the ATTRACT mononucleotide library to build missing atoms for nucleotides", action="store_true")
parser.add_argument("--pdb2pqr",help="Use PDB2PQR to complete missing atoms. If no reference has been specified, analyze the hydrogens to determine histidine/cysteine states", action="store_true")
parser.add_argument("--termini",help="An N-terminus and a C-terminus (5-terminus and 3-terminus for nucleic acids) will be added for each chain", action="store_true")
parser.add_argument("--nter", "--nterm" , dest="nter",
                    help="Add an N-terminus (5-terminus for nucleic acids) for the specified residue number", action="append",
                    type=int, default=[])
parser.add_argument("--cter","--cterm", dest="cter",
                    help="Add a C-terminus (3-terminus for nucleic acids) for the specified residue number", action="append",
                    type=int, default=[])
parser.add_argument("--manual",help="""Enables manual mode.
In automatic mode (default), aareduce tries to produce a PDB file that can be used directly by ATTRACT. In case of missing atoms, a number of
last-resort fixes are attempted that add pseudo-hydrogens at the position of its connected heavy atom. If there are other missing atoms,
an exception is raised.
In manual mode, last-resort fixes are disabled, and missing atoms are simply printed as XXXXXXX in their coordinates. These coordinates
cannot be read by ATTRACT, they need to be edited manually by the user.
""", action="store_true")
parser.add_argument("--trans", "--transfile",dest="transfile",help="Additional trans file that contains additional user-defined atom types (e.g. modified amino acids)", action="append",default=[])
parser.add_argument("--top", "--topfile",dest="topfile",help="Additional topology file in CNS format that contains additional user-defined atom types (e.g. modified amino acids)", action="append",default=[])
parser.add_argument("--patch",dest="patches",
                    help="Provide residue number and patch name to apply", nargs=2, action="append",default=[])
parser.add_argument("--chain", help="Set the chain in the output PDB", default=" ")
parser.add_argument("--startres", help="Set residue number of first residue", type=int, default=1)
parser.add_argument("--startatom", help="Set atom number of first atom", type=int, default=1)
parser.add_argument("--mutate", dest="mutatefiles",
                    help="Provide a 2-column residue mutation file", action="append",default=[])
parser.add_argument("--modres",
                    help="Interpret HETATM records as ATOM if they have a protein backbone", action="store_true")
parser.add_argument("--modbase",
                    help="Interpret HETATM records as ATOM if they have at least three sugar atoms", action="store_true")
parser.add_argument("--batch",
                    help="run aareduce in batch mode. Input and output must be two existing lists of PDBs", action="store_true")
parser.add_argument("--dumppatch",
                    help="Dump all applied patches to a file", action="store_true")
parser.add_argument("--readpatch",
                    help="Read previously applied patches from a file (requires converted input pdb)", action="store_true")
if has_argparse:
    args = parser.parse_args()
else:
    args, positional_args = parser.parse_args()
    args.pdb = None
    args.output = None
    if positional_args:
        args.pdb = positional_args[0]
        if len(positional_args) > 1: args.output = positional_args[1]

assert len(args.chain) == 1, args.chain

if args.rna and args.dna:
    raise ValueError("--dna and --rna are mutually incompatible")

if args.nalib:
    if not args.rna and not args.dna:
        raise ValueError("--nalib requires option --rna or --dna")
    libname = "rnalib"
    if args.dna:
        libname = "dnalib"

if args.heavy and (args.autorefe or args.refe):
    raise ValueError("--(auto)refe and --heavy are mutually incompatible")

if args.autorefe and args.refe:
    raise ValueError("--autorefe and --refe are mutually incompatible")
if args.autorefe:
    args.refe = args.pdb

if args.batch and args.output is None:
    raise ValueError("--batch requires a file list as output argument")

if args.readpatch and len(args.patches):
    raise ValueError("--readpatch and explicit patch specification via --patch are mutually incompatible")

for f in args.topfile:
    assert os.path.exists(f), f
    topfiles.append(f)
for f in args.transfile:
    transfiles.append(f)

topologies = []
for f in topfiles:
    try:
        topologies.append(topology.load(json.load(open(f))))
    except:
        print(f, file=sys.stderr)
        raise
top_residues, top_patches = topology.merge(topologies)

for f in transfiles:
    name = os.path.splitext(os.path.split(f)[1])[0]
    parse_transfile(f, name)

for f in args.mutatefiles:
    assert os.path.exists(f), f
    for l in open(f):
        h = l.find("#")
        if h != -1: l = l[:h]
        ll = l.split()
        if len(ll) == 0: continue
        assert len(ll) == 2, l
        mutations[ll[0]] = ll[1]

if args.nalib:
    nalib = pdbcomplete.load_nalib(libname)

def run(pdbfile):
    pdb = read_pdb(open(pdbfile), pdbfile, add_termini=args.termini,modbase=args.modbase,modres=args.modres)
    pdblines = write_pdb(pdb, args.chain)[0]

    termini_pdb(pdb, args.nter, args.cter)
    patches = {}
    if args.readpatch:
        indata = open(os.path.splitext(pdbfile)[0]+'.patch').readlines()
        indata = [line.split() for line in indata]
        args.patches = indata

    for p in args.patches:
        resid = p[0].strip()
        resindices = [ri for ri,r in enumerate(pdb) if r.resid.strip() == resid]
        if len(resindices) == 0:
            raise ValueError("No residues have resid %s" % resid)
        elif len(resindices) > 1:
            raise ValueError("Multiple residues have resid %s" % resid)
        resid2 = pdb[resindices[0]].resid
        if resid2 not in patches: patches[resid2] = []
        pname = p[1].lower()
        if pname == "none": pname = None
        patches[resid2].append(pname)
    patch_pdb(pdb, patches)

    if args.refe:
        refe = read_pdb(open(args.refe), args.refe, add_termini=args.termini)
        patch_pdb(refe, patches)
        if not args.heavy:
            update_patches(refe, top_patches)
        set_reference(pdb, refe)
    if args.nalib:
        pdbcomplete.apply_nalib(pdb, nalib, args.heavy, args.manual)
    if args.pdb2pqr:
        pdblines = write_pdb(pdb, args.chain, one_letter_na = True)[0]
        pqrlines = run_pdb2pqr(pdblines)
        pqr = read_pdb(pqrlines, "<PDB2PQR output from %s>" % pdbfile)
        pdbcomplete.pdbcomplete(pdb, pqr)
        if not args.heavy and not args.refe:
            update_patches(pdb, top_patches)

    if args.refe:
        pdbfix(pdb, refe)

    if not args.manual:
        pdb_lastresort(pdb)
        check_pdb(pdb, heavy=args.heavy)
    pdblines, mapping = write_pdb(pdb, args.chain, heavy=args.heavy)
    return pdblines, mapping, pdb

if args.batch:
    infiles = read_filelist(args.pdb)
    for f in infiles:
        assert os.path.exists(f), f
    outfiles = read_filelist(args.output)
    for pdb, outfile in zip(infiles, outfiles):
        pdblines, mapping, pdbtop = run(pdb)
        outf = open(outfile, "w")
        for l in pdblines:
            print(l, file=outf)
        outf.close()
        if args.dumppatch:
            outfilep = os.path.splitext(outfile)[0] + ".patch"
            outf = open(outfilep, "w")
            for i,res in enumerate(pdbtop):
                if len(res.topology.patches):
                    for p in res.topology.patches:
                        outf.write(str(i+args.startres)+' '+p+'\n')
            outf.close()
else:
    outfile = os.path.splitext(args.pdb)[0] + "-aa.pdb"
    if args.output is not None:
        outfile = args.output
    pdblines, mapping, pdbtop = run(args.pdb)
    outf = open(outfile, "w")
    for l in pdblines:
        print(l, file=outf)
    for v1, v2 in mapping:
        print(v1, v2)
    outf.close()
    if args.dumppatch:
        outfilep = os.path.splitext(outfile)[0] + ".patch"
        outf = open(outfilep, "w")
        for i,res in enumerate(pdbtop):
            if len(res.topology.patches):
                for p in res.topology.patches:
                    outf.write(str(i+args.startres)+' '+p+'\n')
        outf.close()
