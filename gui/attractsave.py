import os
import Spyder

def make_relpath(outpdir, m):
  try:
    items = [a for a in m]
    items = enumerate(items)
    ar = True
  except TypeError:
    items = list(m.__dict__.items())
    ar = False
  for k,v in items:
    if isinstance(v, Spyder.File): 
      nam = v.name
      rel = os.path.relpath(nam, outpdir)
      if len(os.path.split(rel)[0]) == 0:
        vv = type(v)(
         rel,
         fileformat=v.fileformat(),
         mode=v.mode,
         format=v.format(),
        ) 
        if ar:
          m[k] = vv
        else:
          setattr(m,k,vv) 
    elif isinstance(v, Spyder.Filename):	  
      nam = v.name
      rel = os.path.relpath(nam, outpdir)
      if len(os.path.split(rel)[0]) == 0:
        vv = type(v)(rel)
        if ar:
          m[k] = vv
        else:
          setattr(m,k,vv) 

    elif isinstance(v, Spyder.Object):
      make_relpath(outpdir, v)      
      
      
def save(m, outp, *args):
  v = m.get() 
  outpdir = os.path.split(outp)[0]  
  if v is not None:
    make_relpath(outpdir, v)
    v.tofile(outp)
    print("Form object saved to %s" % outp)
  else:
    print("Cannot save form: does not contain a valid object")

def generate(m, outp, *args):
  v = m.get() 
  outpdir = os.path.split(outp)[0]
  if v is not None:
    make_relpath(outpdir, v)
    sh = os.path.splitext(outp)[0] + ".sh"
    script = v.generate()
    fsh = open(sh, "w")
    fsh.write(script+"\n")
    fsh.close()
    os.system("chmod +x %s" % sh)
    print("Shell script generated: %s" % sh)
  else:
    print("Cannot generate shell script: form does not contain a valid object")
