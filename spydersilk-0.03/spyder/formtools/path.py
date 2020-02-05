import os

def make_relpath(outpdir, m):
  from Spyder import File, Resource, Object, Filename, String, Data
  ar, dic = False, False
  if isinstance(m, dict):
    dic = True
    items = m.items()
  elif isinstance(m, String) or isinstance(m, Data):
    return
  else:
    try:
      items = [a for a in m]
      items = enumerate(items)
      ar = True
    except TypeError:
      items = list(m.__dict__.items())    
  for k,v in items:
    if isinstance(v, Resource): 
      if v.filename is None: continue
      nam = v.filename      
      rel = os.path.relpath(nam, outpdir)
      if rel.startswith(".."): #os.path.relpath does this for /tmp, very annoying
        rel = nam
      v.filename = rel
    elif isinstance(v, File): 
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
        if ar or dic:
          m[k] = vv
        else:
          setattr(m,k,vv) 
    elif isinstance(v, Filename):          
      nam = v
      rel = os.path.relpath(nam, outpdir)
      if rel.startswith(".."): #os.path.relpath does this for /tmp, very annoying
        rel = nam      
      if len(os.path.split(rel)[0]) == 0:
        vv = type(v)(rel)
        if ar or dic:
          m[k] = vv
        else:
          setattr(m,k,vv) 

    elif isinstance(v, Object):
      make_relpath(outpdir, v)      

def make_abspath(m):
  from Spyder import File, Resource, Object, Filename, String, Data
  ar, dic = False, False
  if isinstance(m, dict):
    dic = True
    items = m.items()
  elif isinstance(m, String) or isinstance(m, Data):
    return    
  else:
    try:
      items = [a for a in m]
      items = enumerate(items)
      ar = True
    except TypeError:
      items = list(m.__dict__.items())    
  for k,v in items:
    if isinstance(v, Resource): 
      if v.filename is None: continue
      nam = v.filename      
      absp = os.path.abspath(nam)
      v.filename = absp
    elif isinstance(v, File): 
      nam = v.name
      absp = os.path.abspath(nam)
      vv = type(v)(
        absp,
        fileformat=v.fileformat(),
        mode=v.mode,
        format=v.format(),
      ) 
      if ar or dic:
        m[k] = vv
      else:
        setattr(m,k,vv) 
    elif isinstance(v, Filename):          
      nam = v.name
      absp = os.path.abspath(nam)
      vv = type(v)(absp)
      if ar or dic:
        m[k] = vv
      else:
        setattr(m,k,vv) 

    elif isinstance(v, Object):
      make_abspath(v)      
