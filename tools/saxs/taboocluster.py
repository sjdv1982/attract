import sys, numpy as np, weakref
import itertools
import os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
from _read_struc import read_struc

maxcollectcluster = 10000
MINCHILD = 100
MAXCHUNK = 40000
CLUSTERING = [30, 10, 8, 6, 5, 4, 3.5, 3, 2, 0.1, 0]
#last CLUSTERING must be 0, second-to-last is deredundant criterion
CLUST_MARGIN = 2 #2 is the theoretical worse case; change to 9999 to disable all effects of clustering
#CLUST_MARGIN = 9999

chicutoff = float(sys.argv[1])
chicutoffmin = chicutoff-1.0
max_rmsd = float(sys.argv[2])
maxstruc = int(sys.argv[3])
max_msd = max_rmsd**2
maxclusterlevel = CLUSTERING.index(max_rmsd)
coor_structures = sys.argv[4]
coor_clusters = sys.argv[5]
chi_structures = sys.argv[6]
chi_clusters_filename = sys.argv[7]
dat_structures = sys.argv[8]
dat_clusters_filename = sys.argv[9]
header_structures, dof_structures = read_struc(dat_structures)
header_clusters, dof_clusters = read_struc(dat_clusters_filename)
dof_structures = list(dof_structures)
dof_clusters = list(dof_clusters)
coor_structures = np.load(coor_structures)
coor_clusters = np.load(coor_clusters)
chi_structures = np.load(chi_structures)
chi_clusters = np.load(chi_clusters_filename)
for a in coor_clusters, coor_structures: assert len(a.shape) == 2 and a.shape[1] % 3 == 0, a.shape

arr = [coor_structures, coor_clusters]
nstruc = []
for anr, a in enumerate(arr): 
  ncoor = a.shape[1] / 3
  arr[anr] = a.reshape(a.shape[0], ncoor, 3)
  nstruc.append(len(a))
coor_structures, coor_clusters = arr

ranks = [np.arange(len(s))+1 for s in arr]

   
MAX_SIZEKEY = 100000000000    
MAX_CLUSTERING = len(CLUSTERING) - 1
class Cluster(object):
  __slots__ = ("clustid", "_splittable", "clusterlevel", "coors", "ranks", "all_ranks", "children", "nodes", "totstruc", "parent", "connections", \
    "back_connections", "_splitting", "_checking_delete", "conlevel")
  def __init__(self, clustid, clusterlevel, coors, ranks):
    self.clustid = clustid
    self._splittable = True
   
    if coors is not None and len(coors) == 1: #singleton
      self.clusterlevel = MAX_CLUSTERING      
    else:  
      assert clusterlevel < MAX_CLUSTERING
      self.clusterlevel = clusterlevel #contains the clusterlevel of the cluster itself, not of the cluster children!
    if self.clusterlevel == MAX_CLUSTERING:
      self._splittable = False  
    self.coors = coors
    self.ranks = ranks
    self.all_ranks = set(ranks)
    self.children = []
    self.nodes = 1
    if coors is not None:
      self.totstruc = coors.shape[0]
    self.parent = None
    self.connections = []
    self.back_connections = [] 
    self._splitting = False
    self._checking_delete = False
    self.conlevel = self.clusterlevel
    if self.clusterlevel is None:
      self.conlevel = -1
    if coors is not None and clusterlevel == MAX_CLUSTERING - 1: #deredundant level
      r = self
      for cnr in range(len(self.coors)):
        c = Cluster(self.clustid + (cnr,), MAX_CLUSTERING, coors[cnr:cnr+1], ranks[cnr:cnr+1])
        c.parent = r
        self.children.append(c)
      self.nodes = len(self.children)  
      self._splittable = False

  def _cluster(self, clusterlevel):
    assert not self.children
    
    c = self.coors
    #clus: coordinates of the first structure of each cluster    
    clus = c[:1]
    #clus_indices: the coors indices of the structures of each cluster 
    clus_indices = [[0]]
    chunksize = 20
    
    radius = CLUSTERING[clusterlevel]
    max_sd = radius * radius * c.shape[1]
    
    #This variable keeps track, for every structure in the chunk, into which new cluster it is sorted
    which_new_clust = np.zeros(chunksize, dtype=int)
    
    clustid = self.clustid
    if clustid is None: clustid = ()
    for n in range(1, len(c), chunksize):
      chunk = c[n:n+chunksize]
      
      #intra-chunk msd
      d = chunk[:, np.newaxis, :, :] - chunk[np.newaxis, :, :, :]
      intra_msd = np.einsum("...ijk,...ijk->...i", d,d)
      
      #chunk-cluster msd 
      d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
      inter_msd = np.einsum("...ijk,...ijk->...i", d,d)
      
      """
      tmp = np.sqrt(inter_msd/c.shape[1])
      for i in range(chunksize):
	for j in range(len(clus)):
	  print "RMSD", i+n+1, clus_indices[j][0]+1, tmp[i][j]
      """	  
      
      for nn in range(len(chunk)):
        sort_clust = None        
        close_intra_clusts = (intra_msd[nn] < max_sd).nonzero()[0]        
        intra_new_clusts = [which_new_clust[k] for k in close_intra_clusts if k < nn and which_new_clust[k] != -1]
        if len(intra_new_clusts):          
          sort_clust = min(intra_new_clusts)
        close_inter_clusts = (inter_msd[nn] < max_sd).nonzero()[0]
        if len(close_inter_clusts):
          sort_clust2 = min(close_inter_clusts)
          if sort_clust is None or sort_clust > sort_clust2:
            sort_clust = sort_clust2
        if sort_clust is None:
          #new cluster
          #print "NEWCLUST", nn+n+1
          which_new_clust[nn] = len(clus)
          clus = np.append(clus, chunk[nn][np.newaxis,:,:], axis=0)
          clus_indices.append([n+nn])
        else:
          clus_indices[sort_clust].append(n+nn)
          which_new_clust[nn] = -1
    
    indices = [i[0] for i in clus_indices]
    
    #Re-cluster to the lowest RMSD cluster
    clus_indices = [i[:1] for i in clus_indices]
    for n in range(0, len(c), chunksize):
      chunk = c[n:n+chunksize]      
      d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
      inter_msd = np.einsum("...ijk,...ijk->...i", d,d)      
      sort_clusts = np.argmin(inter_msd, axis=1)       
      for nn in range(len(chunk)):
        if (n+nn) in indices: continue 
        sort_clust = sort_clusts[nn]
        clus_indices[sort_clust].append(n+nn)

    for cnr,c in enumerate(clus_indices):
      ind = clus_indices[cnr]
      coors = self.coors[ind]
      ranks = self.ranks[ind]
      c = Cluster(clustid+(cnr+1,), clusterlevel, coors, ranks)
      self.children.append(c)
    self.nodes = len(self.children)
    self.totstruc = sum([c.totstruc for c in self.children])    
    return clus, indices
  def cluster(self, clusterlevel):
    assert clusterlevel < MAX_CLUSTERING
    clus, indices = self._cluster(clusterlevel)
    self.coors = clus
    self.ranks = self.ranks[indices]
    r = self
    for c in self.children: 
      c.parent = r
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False
  def dissolve(self, clusterlevel):
    self.ranks = np.concatenate([c.ranks for c in self.children])
    newchildren = []
    coors = []
    while len(self.children):
      child = self.children.pop()
      child.cluster(clusterlevel)      
      coors.append(child.coors)
      newchildren += child.children
      self.totstruc = child.totstruc
    self.coors = np.concatenate(coors, axis=0)
    self.children = newchildren
    r = self
    for c in self.children: 
      c.parent = r
    self.nodes = len(newchildren)
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False  
  def reorganize(self):
    for c in self.children:
      c.reorganize()
    if not len(self.children): return
    if len(self.children) >= MINCHILD: return
    oldchildren = [c for c in self.children if len(c.children)]
    if not len(oldchildren): return       
    while len(self.children) < MINCHILD and len(oldchildren):
      child = oldchildren.pop(0)
      self.children.remove(child)
      self.children += child.children
      oldchildren += [c for c in child.children if len(c.children)]    
    r = self
    coors = [c.coors[0] for c in self.children]
    self.coors = np.array(coors)
    for c in self.children: 
      c.parent = r
    self.nodes = sum([c.nodes for c in self.children])
    for c in self.children:
      if c._splittable:
        break
    else:
      self._splittable = False  

  def split(self):
    if not self.children:
      assert self.clusterlevel is not None
      if self.clusterlevel == MAX_CLUSTERING: return False      
      self.clusterlevel += 1      
      if self.clusterlevel == MAX_CLUSTERING: 
        self._splittable = False
        return False
      self.cluster(self.clusterlevel)
      r = self
      while len(self.children) == 1:        
        child = self.children[0]
        self.clusterlevel = child.clusterlevel
        if self.clusterlevel >= MAX_CLUSTERING - 1: 
          self._splittable = False
          break
        self.clusterlevel += 1
        child.cluster(self.clusterlevel)
        children = child.children
        for c in children: c.parent = r
        self.coors = child.coors
        self.ranks = child.ranks
        self.children = children
        self.nodes = len(children)
        self.totstruc = sum([c.totstruc for c in children])
        
      self.parent.add_nodes(self.nodes - 1)
      for c in self.children:
        if c._splittable:
          break
      else:
        self._splittable = False      
      return True
    else:
      self._splitting = True
      oldnodes = self.nodes
      ok = False
      for c in self.children:
        c_splittable = c._splittable
        has_split = c.split()
        if has_split: ok = True
      newnodes = self.nodes
      self._splitting = False
      if self.parent is not None and newnodes > oldnodes:
        self.parent.add_nodes(newnodes-oldnodes)
      for c in self.children:
        if c._splittable:
          break
      else:
        self._splittable = False
      return ok
      
  def add_nodes(self, nodes):    
    self.nodes += nodes
    if self.parent is not None and not self._splitting:
      self.parent.add_nodes(nodes)

  def check_deletion(self):
    if self._checking_delete: return
    self._checking_delete = True
    while 1:
      has_c1, has_c2 = len(self.back_connections), len(self.connections)
      if has_c1 and has_c2: break
      pos = self.clustid[0]
      if pos == 0 and has_c2: break
      if pos == 1 and has_c1: break
      if self not in clusters[pos]: break
      clusters[pos].remove(self)
      for o in list(self.back_connections):
        o.connections.remove(self)
        o.check_deletion()
      for o in list(self.connections):
        o.back_connections.remove(self)
        o.check_deletion()
      break  
    self._checking_delete = False
    
    
  def decompose(self, fwd):
    c1 = self.coors    
    if fwd:
      if not self.connections: return
      others = list(self.connections)
      c2 = np.array([c.coors[0] for c in self.connections])      
    else:
      if not self.back_connections: return
      others = list(self.back_connections)
      c2 = np.array([c.coors[0] for c in self.back_connections])      
    
    chunksize = MAXCHUNK/len(c1)
    for chunkpos in range(0, len(c2), chunksize):
      c2_chunk = c2[chunkpos:chunkpos+chunksize]
      others_chunk = others[chunkpos:chunkpos+chunksize]
      d = c1[:, np.newaxis, :, :] - c2_chunk[np.newaxis, :, :, :]

      o_max_rmsd = []
      for o in others_chunk:
        if o.clusterlevel is None:
          mx = 10**6
        else:
          mx = CLUSTERING[o.clusterlevel] * CLUST_MARGIN
        o_max_rmsd.append(mx)
      c_max_rmsd = []
      for child in self.children:
        mx = CLUSTERING[child.clusterlevel] * CLUST_MARGIN
        c_max_rmsd.append(mx)
      
      max_rmsd0 = np.array(c_max_rmsd)[:, np.newaxis] + np.array(o_max_rmsd)[np.newaxis, :]
      max_rmsd2 = max_rmsd0 + max_rmsd
      max_sd = (max_rmsd2**2) * c1.shape[1]
          
      msd = np.einsum("...ijk,...ijk->...i", d,d)
      if fwd:
        ocon = [o.back_connections for o in others_chunk]
        childcon = [c.connections for c in self.children]
      else:
        ocon = [o.connections for o in others_chunk]
        childcon = [c.back_connections for c in self.children]
      for o in ocon:
        o.remove(self)
      
      for childnr, onr in zip(*np.where(msd < max_sd)):
        #print childnr, onr, len(self.children), len(ocon)
        ocon[onr].append(self.children[childnr])
        childcon[childnr].append(others_chunk[onr])
          
      for o in others:
        o.check_deletion()
      
    
  def verify(self):
    if len(self.children): return
    if self.totstruc > 1: return
    cons = [con for con in self.connections if not len(con.children) and con.totstruc == 1]
    if not len(cons): return
    concoors = np.concatenate([con.coors[:1] for con in cons])
    c2 = self.coors[0]
    c1 = concoors
    d = c1 - c2
    msd = np.einsum("...ijk,ijk->...i", d,d)
    msd_low = (msd<(max_rmsd**2*c2.shape[0]))
    for n in range(len(cons)):
      assert msd_low[n]
    return

  def all_children(self):
    if not len(self.children):
      yield self
    else:  
      for cc in self.children:
        for v in cc.all_children():
          yield v        

  def all_children_of_clusterlevel(self,clusterlevel):
    cl = None
    if self.parent is not None: cl = self.parent.clusterlevel 
    #print "CLUST!", self.minclusterlevel, self.clusterlevel, cl
    if self.clusterlevel is not None and self.clusterlevel >= clusterlevel and (cl is None or cl <= clusterlevel):# 
      yield self
    else:  
      for cc in self.children:
        for v in cc.all_children_of_clusterlevel(clusterlevel):
          yield v        

#Build cluster tree
clusters = []    
for n, atoms in enumerate((coor_structures, coor_clusters)):    
  c = Cluster((n,), None, atoms, ranks[n])
  

  clusterlevel = 0
  c.cluster(clusterlevel)
  for clusterlevel in range(1, maxclusterlevel+1):
    #if you are interested in 6 A clusters, don't dissolve that clusterlevel!
    if len(c.children) >= MINCHILD: break
    c.dissolve(clusterlevel)
  count = 0
  assert c.clusterlevel is None or c.clusterlevel == MAX_CLUSTERING
  
  def split_all(c):
    global count
    if not c._splittable: return
    if not len(c.children):
      ok = c.split()
      count += 1
      if not (count % 500): print >> sys.stderr, n+1, count,  "/", len(atoms[n])
      if not ok: return
    for cc in c.children:
      split_all(cc)
  
  split_all(c) 
  c.reorganize() 
  print >> sys.stderr, n+1, nstruc[n], CLUSTERING[clusterlevel], len(c.children), c.nodes
  assert c.clusterlevel is None or c.clusterlevel == MAX_CLUSTERING
  clusters.append([c])


from copy import deepcopy
structuretree = deepcopy(clusters[0])
sixclust = []
count = 0
for s in structuretree:
  #print s.clusterlevel, len(s.children), [child.clusterlevel for child in s.children]
  for c in s.all_children_of_clusterlevel(maxclusterlevel):
    cl = None
    if c.parent is not None: cl = c.parent.clusterlevel
    count += 1
  sixclust += list(s.all_children_of_clusterlevel(maxclusterlevel))


#Initialize tree connections
c1, c2 = clusters[0][0], clusters[1][0]
for cc in c1.children:
  cc.connections = [c2]
  c2.back_connections.append(cc)
clusters[0] = c1.children  
  

def decompose(clusnr):  
  best_conlevel = None
  best = None
  clus = clusters[clusnr]
  for nn in range(len(clus)):
    c = clus[nn]
    if not len(c.children): continue
    if best_conlevel is None or c.conlevel > best_conlevel:
      best_conlevel = c.conlevel
      best = nn
  if best is None: return False
  c = clus[best]
  clusters[clusnr].remove(c)
  if (clusnr == 0):
    c.decompose(fwd=True)
  else:
    c.decompose(fwd=False)
  if not (count % 100): print >> sys.stderr, clusnr+1, count,  "/", len(arr[clusnr])
  for child in c.children:
    conlevel = 0
    for con in itertools.chain(child.back_connections,child.connections):
      v = con.clusterlevel
      if v is not None and v > conlevel: conlevel = v
    child.conlevel = conlevel * child.clusterlevel
  clusters[clusnr] += c.children
  
  for cc in c.children:
    cc.check_deletion()
  return True

count = 0
done1, done2 = False, False
while not (done1 and done2):
  count += 1
  if not done1: 
    done1 = not decompose(0)
  if not done2: 
    done2 = not decompose(1)
    
        
print >> sys.stderr, [len(c) for c in clusters]
print >> sys.stderr, [sum([cc.nodes for cc in c]) for c in clusters]

#Verification
for c in clusters:
  for cc in c:
    cc.verify()

inter = []
for c1 in clusters[0]:
  r1 = c1.ranks[0]
  for c2 in c1.connections:
    r2 = c2.ranks[0]
    inter.append((r1,r2))
    
swapped = set()
#check if new structures could replace old cluster heads
for r1, r2 in inter:
  if chi_structures[r1-1] < chi_clusters[r2-1] and r1 not in swapped:
    swapped.add(r1)#swap at most one time
    chi_clusters[r2-1] = chi_structures[r1-1]
    dof_clusters[r2-1] = dof_structures[r1-1]
    
    
inter = set([r1 for r1,r2 in inter])
#add new clusters and remove the structures which are similar to them
for c in sixclust:
  tmp = c.ranks[0]
  tmpchi = chi_structures[tmp-1]
  if tmp not in inter and tmpchi < chicutoff:
    #add as cluster head
    worst = np.amax(chi_clusters)
    added = False
    if len(chi_clusters) < maxcollectcluster:
      chi_clusters = np.append(chi_clusters,tmpchi)
      dof_clusters.append(dof_structures[tmp-1])
      added = True
      
    elif tmpchi < worst and worst > chicutoffmin:
      replace = np.argmax(chi_clusters)
      chi_clusters[replace] = tmpchi
      dof_clusters[replace] = dof_structures[tmp-1]
      added = True
      
    elif worst < chicutoffmin and tmpchi < chicutoffmin:
      chi_clusters = np.append(chi_clusters,tmpchi)
      dof_clusters.append(dof_structures[tmp-1])
      added = True
      
    if added:
      for j in c.ranks:
      #remove
	inter.add(j)



#Save all data
out = open(dat_clusters_filename,'w')  

for hh in header_clusters:
  out.write(hh+'\n')
  
for nstruc,ll in enumerate(dof_clusters):
  out.write('#'+str(nstruc+1)+'\n')
  l1, l2 = ll
  for l in l1:
    out.write(l+'\n')
    
  for l in l2:
    out.write(l+'\n')
    
out.close()
np.save(chi_clusters_filename,chi_clusters)
for hh in header_structures:
  print hh
  
nstruc = 0
for j,ll in enumerate(dof_structures):
  if j+1 in inter: continue
  nstruc+=1
  print "#"+str(nstruc)
  l1, l2 = ll
  for l in l1:
    print l
    
  for l in l2:
    print l

  
    


  


