import numpy as np
def read_situs(situsfile):
  header = open(situsfile).readline()
  h = header.split()
  assert len(h) == 7, header
  voxelsize = float(h[0])
  origin = np.array([float(v) for v in h[1:4]])
  dimensions = tuple([int(v) for v in h[4:7]])
  data = np.genfromtxt(situsfile,dtype="float32",skip_header=2,skip_footer=1)
  nvox = dimensions[0]*dimensions[1]*dimensions[2]
  for lastline in open(situsfile).xreadlines():
    pass
  lastdata = np.array([float(v) for v in lastline.split()])
  data = data.reshape(-1)
  data = np.append(data, lastdata)
  data = data.reshape(dimensions, order='F')
  data = np.ascontiguousarray(data)
  return data, voxelsize, origin

def write_situs(situsfile, data, voxelsize, origin):
  f = open(situsfile, "w")
  f.write("%.6f " % voxelsize)
  f.write("%.6f %.6f %.6f " % tuple(origin))
  f.write("%d %d %d\n\n" % data.shape)
  d = data.reshape(-1, order='F')
  count = 0
  for v in d:
    f.write("%11.6f " % v)
    count += 1
    if not count % 10: f.write("\n")
    #else: f.write(" ")
  f.close()
