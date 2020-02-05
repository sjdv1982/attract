def check_embedded(m, attrname = ""):
  import Spyder
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

def embed(m):
  import Spyder
  try:
    items = [a for a in m]
    items = enumerate(items)
    ar = True
  except TypeError:
    items = list(m.__dict__.items())
    ar = False
  for k,v in items:
    if isinstance(v, Spyder.Resource): 
      v.embed()
    elif isinstance(v, Spyder.Object):
     embed(v)      
      