from .model import model
from .controller import controller
from .generate_typetree import generate_typetree, _generate_typetree
from .path import make_abspath, make_relpath
from .embed import embed, check_embedded
from .controller import controller_leaf

def apply_delta(delta, controller, mode=None):
  import spyder, Spyder
  if delta is None: return
  is_resource = controller._model()._is_resource
  is_file = controller._model()._is_file
  is_leaf = isinstance(controller, controller_leaf)  
  if is_resource:
    if mode == "json" and not isinstance(delta, str) and not isinstance(delta, dict): 
      delta = tuple(delta)
    elif isinstance(delta, dict):   
      delta = controller._model()._resourcetype.fromdict(delta)          
    else:
      delta = controller._model()._resourcetype(delta)
    controller._set(delta)
    return
  if is_leaf:
    if is_file and mode is None and isinstance(delta, dict):   
      delta = Spyder.File(delta)
    controller._set(delta)
    return
  
  if isinstance(delta, dict):
    for k in delta:
      subdelta = delta[k]
      subcontroller = getattr(controller, k)
      apply_delta(subdelta, subcontroller, mode)
  elif isinstance(delta, list):
    for n,subdelta in enumerate(delta):
      subcontroller = controller[n]
      apply_delta(subdelta, subcontroller, mode)
  else:
    controller._set(delta)

from . import cgimodule as cgi
