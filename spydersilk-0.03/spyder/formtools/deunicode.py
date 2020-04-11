import sys
python3 = (sys.version_info[0] == 3)

def _deunicode_dict(dic):  
  ret = {}
  for k,v in dic.items():
    ret[k.encode()] = _deunicode(v)
  return ret

def _deunicode_list(l):  
  return [_deunicode(ll) for ll in l]
  
def _deunicode(obj):
  if isinstance(obj, dict): return _deunicode_dict(obj)
  elif isinstance(obj, list): return _deunicode_list(obj)
  elif isinstance(obj, unicode): return obj.encode()
  else: return obj
  
def deunicode(obj):
  if python3: return obj 
  return _deunicode(obj)
  