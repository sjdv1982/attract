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
    if isinstance(v, Spyder.Resource): 
      if v.filename is None: continue
      nam = v.filename      
      rel = os.path.relpath(nam, outpdir)
      if rel.startswith(".."): #os.path.relpath does this for /tmp, very annoying
        rel = nam
      v.filename = rel
    elif isinstance(v, Spyder.File): 
      nam = v.name
      rel = os.path.relpath(nam, outpdir)
      if rel.startswith(".."): #os.path.relpath does this for /tmp, very annoying
        rel = nam      
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
      if rel.startswith(".."): #os.path.relpath does this for /tmp, very annoying
        rel = nam      
      if len(os.path.split(rel)[0]) == 0:
        vv = type(v)(rel)
        if ar:
          m[k] = vv
        else:
          setattr(m,k,vv) 

    elif isinstance(v, Spyder.Object):
      make_relpath(outpdir, v)      
      
def check_embedded(m, attrname = ""):
  try:
    items = [a for a in m]
    items = enumerate(items)
    ar = True
  except TypeError:
    items = list(m.__dict__.items())
    ar = False
  for k,v in items:
    attrname2 = attrname
    if ar: attrname2 += "[%s]" % k
    elif attrname == "": attrname2 = str(k)
    else: attrname2 += ".%s" % k    
    if isinstance(v, Spyder.Resource): 
      if v._data is not None:
        return attrname2
    elif isinstance(v, Spyder.Object):
      embedded = check_embedded(v, attrname2)      
      if embedded is not None: return embedded
  
def save(m, outp, *args):
  v = m._get()   
  outpdir = os.path.split(outp)[0]  
  if v is not None:
    make_relpath(outpdir, v)
    v.tofile(outp)
    print("Form object saved to %s" % outp)
  else:
    print("Cannot save form: does not contain a valid object")

def generate(m, outp, *args):
  v = m._get() 
  outpdir = os.path.split(outp)[0]
  if v is not None:
    embedded = check_embedded(v)
    if embedded is not None:
      print("Cannot generate shell script: %s is an embedded resource, not a file name" % embedded)
      return
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
