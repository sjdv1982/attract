import sys, os
import spyder, Spyder
import attractmodel
import resources

inputfile = sys.argv[1]
inputdata = open(inputfile).read()
inputdir = os.path.split(inputfile)[0]
if len(inputdir): os.chdir(inputdir)
spydertypename, spyderdict = spyder.core.parse(inputdata)
spydertype = getattr(Spyder, spydertypename)
model = spydertype.fromdict(spyderdict)
resources.embed(model)
print repr(model)