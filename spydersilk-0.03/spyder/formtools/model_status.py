import traceback
import spyder

_undef = "Undefined"
_async_inc = "Out of sync (incomplete)"
_inc = "Incomplete"
_opt_inc = "Incomplete, optional"
_async_inv = "Out of sync (invalid)"
_inv = "Invalid"
_opt_inv = "Invalid, optional"
_async_ech = "Out of sync (in children)"
_ech = "Error in children"
_opt_ech = "Optional; error in children"
_def = "Defined"
_ok = "OK"

def _get_member(member,controller):
  m = member
  try:
    subcontroller = controller._children[member]
    name = getattr(subcontroller._form, "name")
    if name: m = "%s (%s) " % (m, name)
  except:
    pass
  return "Member %s:" % m

def _fully_ok(model, controller):  
  c = model._constructionstate
  opt = model._is_optional
  ok = False
  if c == "full" and not len(model._async_children): ok = True
  if c == "empty" and opt: ok = True
  if not ok: 
    return False
  
  if not hasattr(controller, "_difchildren"): return True
  for childname in controller._difchildren:
    child = model._children[childname]
    controllerchild = controller._children[childname]
    if not _fully_ok(child, controllerchild): 
      return False
  return True  
    

def _formdefined(model, controller):  
  if controller is None: return []
  if not hasattr(controller, "_difchildren"): return []    
  formdefined = []
  for childname in controller._difchildren:
    child = model._children[childname]
    controllerchild = controller._children[childname]
    if not _fully_ok(child, controllerchild): 
      formdefined.append(childname)
  return sorted(formdefined)

def _add_formdefined(model,controller=None):
  ret = []
  formdefined = []
  if controller is not None:
    formdefined = _formdefined(model, controller)
  for m in formdefined:
    child = model._children[m]
    ret.append(_get_member(m, controller))
    childcontroller = None
    if controller is not None:
      childcontroller = controller._children.get(m, None)      
    mstat = _model_status(child,childcontroller)
    for l in mstat:
      ret.append("  " + l)
  return ret
  
def _missing(model):
  if not len(model._requiredmembers): return []
  defined = set()
  if model._model__dif is not None: defined = model._model__dif.keys()
  missing = model._requiredmembers.difference(defined)
  return sorted(list(missing))  

def _add_missing(model,controller=None):
  ret = []
  missing = _missing(model)
  if missing: 
    ret.append("Missing members: %s" % str(missing))
  if controller is not None:
    formdefined = _formdefined(model, controller)
    missing = sorted(list(set(missing+formdefined)))
  for m in missing:
    child = model._children[m]
    childcontroller = None
    if controller is not None:
      childcontroller = controller._children.get(m, None)
    ret.append(_get_member(m, controller))
    mstat = _model_status(child,childcontroller)
    for l in mstat:
      ret.append("  " + l)
  return ret

def _invalid(model):
  invalid = []
  for childname in model._children:
    child = model._children[childname]
    cc = child._constructionstate
    copt = child._is_optional
    if cc == "full" and not len(child._async_children): continue
    if cc == "empty" and copt: continue
    invalid.append(childname)
  return sorted(invalid)
  
def _format_exception(exc):
  if issubclass(exc[0], spyder.ValidationError) or issubclass(exc[0], spyder.ConstructionError):      
    return traceback.format_exception_only(*exc[:2])
  return traceback.format_exception(*exc)
  
def _add_invalid(model,controller=None):
  ret = []
  exc = model._model__exc
  if exc is not None:
    ret.append("Exception:")    
    exc2 = _format_exception(exc)
    ret.append("  " + "\n  ".join(exc2))  
  #if issubclass(exc[0], spyder.ValidationError): 
  #  return ret
  
  invalid = _invalid(model)
  if invalid: 
    ret.append("Invalid members: %s" % str(invalid))
  if controller is not None:
    formdefined = _formdefined(model, controller)
    invalid = sorted(list(set(invalid+formdefined)))
  for m in invalid:
    child = model._children[m]
    ret.append(_get_member(m, controller))
    childcontroller = None
    if controller is not None:
      childcontroller = controller._children.get(m, None)      
    mstat = _model_status(child,childcontroller)
    for l in mstat:
      ret.append("  " + l)
  return ret
  
def _model_status(model, controller=None,toplevelcast=False):
  c = model._constructionstate
  opt = model._is_optional
  valued = model._is_valued
  if c == "full" and len(model._async_children): c = "childerror"
  
  if c == "empty":
    if opt:
      assert not valued, model._path
      return [_undef] + _add_formdefined(model,controller)
    else:
      if valued:
        return [_async_inc] + _add_missing(model,controller)
      else:
        return [_inc] + _add_missing(model,controller)
  elif c == "partial":    
    if valued:
      return [_async_inc] + _add_missing(model,controller)  
    else:
      if opt:
        return [_opt_inc] + _add_missing(model,controller)
      else:
        return [_inc] + _add_missing(model,controller)
  elif c == "failed":
    if valued:      
      return [_async_inv] + _add_invalid(model,controller)
    else:
      if opt:
        return [_opt_inv] + _add_invalid(model,controller)
      else:
        return [_inv] + _add_invalid(model,controller) 
  elif c == "childerror":
    if valued:      
      return [_async_ech] + _add_invalid(model,controller)
    else:
      if opt:
        return [_opt_ech] + _add_invalid(model,controller)
      else:
        return [_ech] + _add_invalid(model,controller)         
  elif c == "full":
    assert valued or model._basic or model._complete(), model._path
    if not toplevelcast:
      return _add_formdefined(model,controller)
    if opt:
      return [_def] + _add_formdefined(model,controller)
    else:
      return [_ok] + _add_formdefined(model,controller)
  else:
    raise ValueError(model._constructionstate)
      
def model_status(model,controller=None): 
  ret = _model_status(model,controller,toplevelcast=True)
  return "Status: " + "\n  ".join(ret)