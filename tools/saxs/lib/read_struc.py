import sys, os, numpy as np
sys.path.append(os.environ["ATTRACTTOOLS"])
from _read_struc import read_struc as _read_struc
from .euler2rotmat import euler2rotmat
from topview import TVArray

def read_struc(datfile, coms):
  header, strucs = _read_struc(datfile)
  pivots0 = []
  eulers0 = []
  positions0 = []
  pivot_auto = None
  centered_receptor, centered_ligands = None, None
  for h in header:
    if h.startswith("#pivot "):
      hh = h.split()
      if hh[1] != "auto":
        assert not pivot_auto
        p = float(hh[2]), float(hh[3]), float(hh[4])
        pivots0.append(p) 
        pivot_auto = False
      else:
        pivot_auto = True
    elif h.startswith("#centered "):  
      hh = h.split()
      if hh[2] == "true":
        v = True
      elif hh[2] == "false":  
        v = False
      else:
        raise ValueError(h)
      assert hh[1] in ("receptor:", "ligands:"), h
      if hh[1] == "receptor:":
        centered_receptor = v
      else:
        centered_ligands = v
  assert centered_receptor is not None
  assert centered_ligands is not None
  assert pivot_auto is not None
  if pivot_auto:
    pivots = np.array(coms, dtype="float32")
  else:    
    pivots = np.array(pivots0, dtype="float32")
  nbodies = len(pivots)
  for s1, s2 in strucs:
    curr_eulers = np.ndarray((nbodies,3), dtype="float32")
    curr_positions = np.ndarray((nbodies,3), dtype="float32")    
    assert len(s2) == nbodies, (len(s2), nbodies)
    for body, l in enumerate(s2):
      dofs = [float(v) for v in l.split()]
      assert len(dofs) == 6 #no modes or ensembles for now
      curr_eulers[body] = dofs[:3]
      curr_positions[body] = dofs[3:]
    eulers0.append(curr_eulers)
    positions0.append(curr_positions)
  nstruc = len(eulers0)
  eulers = np.empty((nstruc,nbodies,3),dtype="float32")
  positions = np.empty((nstruc,nbodies,3),dtype="float32")
  for n in range(nstruc):
    eulers[n] = eulers0[n]
    positions[n] = positions0[n]
    
  return pivots, eulers, positions  
      
def depivotize(pivots, eulers, positions):
  eul = np.ascontiguousarray(np.transpose(eulers, (1,0,2)))
  ret = positions.copy()
  for n in range(len(eul)):
    g_e = TVArray("g_e", eul[n], gpu=True)
    g_rots = TVArray("rot", gpu=True)
    euler2rotmat(g_e) > g_rots
    rots = g_rots.get_data().get()
    rots = rots.reshape(len(rots),3,3)
    p = pivots[n]
    pp = (-p * rots).sum(axis=2) + p
    ret[:,n] += pp
  return ret
