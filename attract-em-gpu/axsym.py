import sys, argparse
import numpy as np

from lib.read_struc import read_struc, depivotize

import _axsym
from math import *



class AxSymAction(argparse.Action):
  def __call__(self, parser, namespace, values, option_string=None):
    if getattr(namespace, self.dest, None) is None:
      setattr(namespace, self.dest, [])
    assert len(values) == 8, values  
    axsym = AxSymmetry(values[0], values[1], 0, *values[2:8])
    assert axsym.symtype >= 0
    getattr(namespace, self.dest).append(axsym)
    
class AxSymmetry(object):
  def __init__(self, ligand, symtype, angle, ax, ay, az, ox, oy, oz):
    self.ligand = int(ligand)
    self.symtype = int(symtype)
    self.angle = float(angle)
    assert symtype > 0
    self.axis = np.array((ax,ay,az),dtype="float")
    self.origin = np.array((ox,oy,oz),dtype="float")


def prepare_symtrans(axsyms, nlig):
  symcopies = np.ones(nlig, "int")      
  for ax in axsyms:
    copies = ax.symtype
    l = ax.ligand - 1;
    if ax.symtype == 0: copies = 2
    symcopies[l] *= copies
  nsymtrans = sum(symcopies) - nlig
  symtrans = _axsym.ffi.new("SymTrans[]", nsymtrans)   

  symcopies = [[n] for n in range(nlig)]     
  pos = 0
  for ax in axsyms:
    copies = ax.symtype
    l = ax.ligand - 1;
    if ax.symtype == 0: copies = 2 
    lc = len(symcopies[l])
    for n in range(1, copies):
      for nn in range(lc):
        symtr = symtrans[pos]
        symtr.ligand = symcopies[l][nn]
        symtr.targetligand = nlig + pos
        symcopies[l].append(symtr.targetligand)
        symtr.origin = list(ax.origin)
        
        if ax.symtype == 0:
          angle = sym.angle / 180 * pi
        else:
          angle = 2.0 * pi / copies * n
      
        #Compute a symmetry matrix from the symmetry axis and the current angle
        rotmatsym = symtr.rotmatsym
        c = cos(angle)
        s = sin(angle)
        t = 1 - c
        x = ax.axis[0]
        y = ax.axis[1]
        z = ax.axis[2]
        rotmatsym[0] = t*x*x + c
        rotmatsym[1] = t*x*y + z*s
        rotmatsym[2] = t*x*z - y*s
        rotmatsym[3] = t*x*y - z*s
        rotmatsym[4] = t*y*y + c
        rotmatsym[5] = t*y*z + x*s
        rotmatsym[6] = t*x*z + y*s
        rotmatsym[7] = t*y*z - x*s
        rotmatsym[8] = t*z*z + c
      
        pos += 1  
  return symtrans   
  

def write_struc(f, pivots, eulers, positions):
  for n in range(len(pivots)):
    print >>f,"#pivot %d %.3f %.3f %.3f" % (n+1, pivots[n][0], pivots[n][1], pivots[n][2])
  print >>f,"#centered receptor: false"
  print >>f,"#centered ligands: false"
  strucnr = 0
  for eu, pos in zip(eulers, positions):
    strucnr += 1
    print >>f,"#" + str(strucnr)
    for eu2, pos2 in zip(eu, pos):
      print >>f," ",  eu2[0], eu2[1], eu2[2], pos2[0], pos2[1], pos2[2]

def apply_axsym(symtrans, eulers, positions):
  nstruc, nlig = eulers.shape[:2]
  arr = np.zeros((6, nstruc, nlig + len(symtrans)))
  args = []
  for n in range(3):
    arr[n,:,:nlig] = eulers[:,:,n]
    args.append(_axsym.ffi.cast("double *", arr[n].ctypes.data))
  for n in range(3):
    arr[3+n,:,:nlig] = positions[:,:,n]
    args.append(_axsym.ffi.cast("double *", arr[3+n].ctypes.data))
  _axsym.lib.apply_axsym(len(symtrans), symtrans, nstruc, nlig + len(symtrans), *args)
  new_eulers = np.transpose(arr[:3], (1,2,0))
  new_positions = np.transpose(arr[3:], (1,2,0))
  return new_eulers, new_positions
  

if __name__ == "__main__":
  a = argparse.ArgumentParser(prog="neo-attract-em.py")
  a.add_argument("datfile")
  a.add_argument("--axsym", nargs="*", default=[], action=AxSymAction)
  args = a.parse_args()

  pivots, eulers, positions0 = read_struc(args.datfile, None) 
  positions = depivotize(pivots, eulers, positions0)  
  
  symtrans = prepare_symtrans(args.axsym, eulers.shape[1])
  new_eulers, new_positions = apply_axsym(symtrans, eulers, positions)
  new_pivots = [(0,0,0)] * len(new_eulers[1])
  write_struc(sys.stdout, new_pivots, new_eulers, new_positions) 
    