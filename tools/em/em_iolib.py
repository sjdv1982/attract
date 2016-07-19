import numpy, struct

def read_situs(situsfile):
  header = open(situsfile).readline()
  h = header.split()
  assert len(h) == 7, header
  voxelsize = float(h[0])
  origin = numpy.array([float(v) for v in h[1:4]])
  dimensions = tuple([int(v) for v in h[4:7]])
  data = numpy.genfromtxt(situsfile,skip_header=1,skip_footer=1)
  nvox = dimensions[0]*dimensions[1]*dimensions[2]
  for lastline in open(situsfile).xreadlines():
    pass
  lastdata = numpy.array([float(v) for v in lastline.split()])
  data = data.reshape(-1)
  data = numpy.append(data, lastdata)
  assert data.size == nvox, (data.size, nvox)
  data = data.reshape(dimensions, order='F')
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

def write_mask(maskfile, data, voxelsize, origin):
  f = open(maskfile, "wb") 
  f.write("ATTRACTMASK")  
  voxelsize = struct.pack("f", voxelsize)
  f.write(voxelsize)  
  origin = struct.pack("fff", *tuple(origin))
  f.write(origin)  
  shape = struct.pack("iii", *data.shape)
  f.write(shape)
  f.write(data.data)
