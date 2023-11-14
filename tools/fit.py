#!/usr/bin/env python3

from __future__ import print_function
import sys

import numpy
from math import *
import rmsdlib

def read_pdb(pdb):
    atoms = []
    lines0 = open(pdb).readlines()
    lines = []
    extralines = []
    for l in lines0:
        if not l.startswith("ATOM"):
            extralines.append((len(lines), l))
            continue
        x = float(l[30:38])
        y = float(l[38:46])
        z = float(l[46:54])
        atoms.append((x,y,z))
        lines.append(l)
    return lines, numpy.array(atoms), extralines

def apply_matrix(atoms, pivot, rotmat, trans):
    ret = []
    for atom in atoms:
        a = atom-pivot
        atom2 = a.dot(rotmat) + pivot + trans
        ret.append(atom2)
    return ret

def cyclesfit(atoms1, atoms2,iterations=5,cutoff=2):
    """
    Performs a PyMol-style fit with rejection cycles
     After every cycle, atoms with more than <cutoff> standard deviations error (in A2)
     get rejected from the selection
    """
    a1, a2 = numpy.array(atoms1), numpy.array(atoms2)
    pivot = numpy.sum(a2,axis=0) / float(len(a2))
    for n in range(iterations):
        rotmat, offset, rmsd = rmsdlib.fit(a1,a2)
        aa2 = apply_matrix(a2, pivot, rotmat, offset)
        dif = aa2 - a1
        d = dif.sum(axis=1)
        sd = d * d
        msd = sd.sum()/len(sd)
        #rmsd = sqrt(msd)
        std = numpy.std(sd)
        if (std < 0.1): std = 0.1

        keep = numpy.less(sd,cutoff*std)
        aa1 = a1[keep]
        aa2 = a2[keep]

        a1 = aa1
        a2 = aa2
    return rmsdlib.fit(a1, a2)

def select_bb(lines, atoms):
    ret = []
    for l, a in zip(lines, atoms):
        if l[13:15] in ("CA","C ","O ","N "): ret.append(a)
    return ret

def select_ca(lines, atoms):
    ret = []
    for l, a in zip(lines, atoms):
        if l[13:15] in ("CA",): ret.append(a)
    return ret

def select_p(lines, atoms):
    ret = []
    for l, a in zip(lines, atoms):
        if l[13:15] in ("P ",): ret.append(a)
    return ret

def select_sel(lines, atoms, selection):
    ret = [ atoms[int(i)-1] for i in selection]
    return ret

def write_pdb(lines, atoms, extralines):
    count = 0
    pos = 0
    while 1:
        while pos < len(extralines):
            p,d = extralines[pos]
            if count < p: break
            print(d.rstrip("\n"))
            pos += 1
        if count == len(lines): break
        l, a = lines[count], atoms[count]
        ll = l[:30] + "%8.3f%8.3f%8.3f" % (a[0],a[1],a[2]) + l[54:].rstrip("\n")
        print(ll)
        count += 1

import sys
import argparse
a = argparse.ArgumentParser(prog="fit.py")
a.add_argument("reference")
a.add_argument("mobile")
a.add_argument("--allatoms", action="store_true")
a.add_argument("--ca", action="store_true")
a.add_argument("--p", action="store_true")
a.add_argument("--rmsd", action="store_true")
a.add_argument("--nofit", action="store_true", help="Calculate RMSD without fitting")
a.add_argument("--iterative", action="store_true")
a.add_argument("--iterative_cycles",type=int,default=5)
a.add_argument("--iterative_cutoff",type=float,default=2)
a.add_argument("--selection1",help="reference atom indices to use for fitting", nargs="+")
a.add_argument("--selection2",help="mobile atom indices to use for fitting", nargs="+")

args = a.parse_args()

#read atoms
lines1, atoms1, extralines1 = read_pdb(args.reference)
lines2, atoms2, extralines2 = read_pdb(args.mobile)

#select backbone
if args.allatoms:
    assert not args.ca
    assert not args.p
    atoms1_fit = atoms1
    atoms2_fit = atoms2
elif args.ca:
    assert not args.p
    atoms1_fit = select_ca(lines1, atoms1)
    atoms2_fit = select_ca(lines2, atoms2)
elif args.p:
    atoms1_fit = select_p(lines1, atoms1)
    atoms2_fit = select_p(lines2, atoms2)
elif args.selection1:
    atoms1_fit = select_sel(lines1, atoms1, args.selection1)
    atoms2_fit = select_sel(lines2, atoms2, args.selection2)
else:
    atoms1_fit = select_bb(lines1, atoms1)
    atoms2_fit = select_bb(lines2, atoms2)
assert len(atoms1_fit) and len(atoms1_fit) == len(atoms2_fit), (len(atoms1_fit), len(atoms2_fit))


if args.nofit:
    d = atoms1_fit - atoms2_fit
    rmsd = numpy.sqrt((d * d).sum() / len(d))
else:
    if args.iterative:
        #perform a Pymol-style iterative fit
        rotmat, offset, rmsd = cyclesfit(atoms1_fit,atoms2_fit, args.iterative_cycles, args.iterative_cutoff)
    else:
        #perform a direct fit
        rotmat, offset, rmsd = rmsdlib.fit(atoms1_fit,atoms2_fit)

pivot = numpy.sum(atoms2_fit,axis=0) / float(len(atoms2_fit))
if args.rmsd or args.nofit:
    print("%.3f" % rmsd, file=sys.stderr)
else:
    fitted_atoms = apply_matrix(atoms2, pivot, rotmat, offset)
    write_pdb(lines2, fitted_atoms, extralines2)
