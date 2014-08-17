import sys, numpy
from em_iolib import read_situs, write_situs, write_mask

situsfile = sys.argv[1]
data, gridspacing, origin = read_situs(situsfile)
print "%.3f" % numpy.sum(data, axis=None)