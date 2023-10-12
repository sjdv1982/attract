import os
import argparse
parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("pdb",help="PDB file to reduce")
parser.add_argument("output",help="all-atom reduced output PDB file", nargs="?")

parser.add_argument("--heavy",help="Ignore all hydrogens", action="store_true")
parser.add_argument("--refe", "--reference",help="Analyze the hydrogens of a reference file to determine histidine/cysteine states")
parser.add_argument("--autorefe",help="Analyze the hydrogens of the input PDB to determine histidine/cysteine states", action="store_true")
parser.add_argument("--dna",help="Automatically interpret nucleic acids as DNA", action="store_true")
parser.add_argument("--rna",help="Automatically interpret nucleic acids as RNA", action="store_true")
parser.add_argument("--nalib",help="Use the ATTRACT mononucleotide library to build missing atoms for nucleotides", action="store_true")
parser.add_argument("--pdb2pqr",help="Use PDB2PQR to complete missing atoms. If no reference has been specified, analyze the hydrogens to determine histidine/cysteine states", action="store_true")
parser.add_argument("--termini",help="An N-terminus and a C-terminus (5-terminus and 3-terminus for nucleic acids) will be added for each chain", action="store_true")
parser.add_argument("--nter", "--nterm" , dest="nter",
                    help="Add an N-terminal H+ (5-terminal phosphate OH for nucleic acids) for the specified residue number", action="append",
                    type=int, default=[])
parser.add_argument("--cter","--cterm", dest="cter",
                    help="Add a C-terminal OH- (3-terminal H for nucleic acids) for the specified residue number", action="append",
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
args = parser.parse_args()

def read_filelist(filelist):
    ret = []
    for l in open(filelist):
        l = l.strip()
        if not len(l): continue
        assert len(l.split()) == 1, (filelist, l)
        ret.append(l)
    return ret

if_results = []
if_files = []
if args.batch:
    infiles = read_filelist(args.pdb)
    outfiles = read_filelist(args.output)
    outpatchfiles = []
    for pdb, outfile in zip(infiles, outfiles):
        if args.dumppatch:
            outfilep = os.path.splitext(outfile)[0] + ".patch"
            outpatchfiles.append(outfilep)
    if_files += infiles
    if_results += outfiles
    if_results += outpatchfiles
else:
    if_files.append(args.pdb)
    outfile = os.path.splitext(args.pdb)[0] + "-aa.pdb"
    if args.output is not None:
        outfile = args.output
    if_results.append(outfile)
    if args.dumppatch:
        outfilep = os.path.splitext(outfile)[0] + ".patch"
        if_results.append(outfilep)

interface = {
    "files": if_files,
    "results": if_results
}

import json
print(json.dumps(interface, indent=2))
