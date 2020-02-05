from functools import partial
from .model import dict2list
from .comparator import comparator
from ..spydercompile import validvar2
from . import arraymanager_dynamic

#"extern" parameter: triggered by a comparator
    

def _get_hardform(controller):
  hardform = controller._hardform
  if hardform == True: return True
  if hardform == False: return False
  assert not controller._toplevel
  return _get_hardform(controller._parent)

def _max_array_child(cont):
  for maxnr in range(len(cont._children),0,-1):
    child = cont._children[maxnr-1]
    if child._always_active: break
    if isinstance(child, controller_leaf) and child._dif(): break
    if isinstance(child, controller) and len(child._difchildren): break
  return maxnr

def _nonevalue(v):
  if v is None or v in (0, 0.0, "", False): return True
  if isinstance(v, tuple) and len(v) == 2 and (v[1] is None or v[1] == "") and validvar2(v[2]):
    return True
  return False

def copy_controller(f1, f2):
  #returns a set of actions
  #performing the actions copies controller f1 onto f2
  
  actions = []
  f2._active = f1._active
  
  if isinstance(f1, controller_leaf):
    f2._value = f1._value
    f2._value2 = f1._value2
    f2._controller_leaf__dif = f1._controller_leaf__dif
        
    c1 = f1._comparator
    c2 = f2._comparator
    c2._value = c1._value
    c2._defined_value = c1._defined_value
        
    assert len(f1._views) == len(f2._views)
    for v1, v2 in zip(f1._views.values(), f2._views.values()):
      actions.append((v2.set, v1.value))
  else:
    f2._difchildren = f1._difchildren
    for child in f1._children:
      ff1 = f1._children[child]
      ff2 = f2._children[child]
      actions += copy_controller(ff1, ff2)
  return actions   
      
def create_controller(
   form, model, typedefault=None,default=None, 
   prename=None,typetree=None, parent=None, membername=None, model_is_func = False
):
    import Spyder, spyder
    assert isinstance(form, spyder.core.spyderform)
    if hasattr(form, "type") and form.type in (None, "none"):
      assert parent is not None
      return  
    if form is None: 
      assert parent is not None
      return
    
    leaf = False
    membernames = form.get_membernames()
    if membernames is None and form.arraycount == 0: 
      leaf = True
    elif hasattr(form, "type"):
      leaf = True
      
    if leaf:
      return controller_leaf(
       form,model,typedefault,default,
       prename,typetree,parent,membername,model_is_func
      )      
    else:
      return controller(
       form,model,typedefault,default,
       prename,typetree,parent,membername,model_is_func
      )    
    
class controller_base(object):
  _listening = False
  def __init__(self, 
   form, model, typedefault=None,default=None,
   prename=None,typetree=None, parent=None, membername=None, model_is_func = False
  ):    
    self._form = form
    if model_is_func:
      self._model = model
    else:
      self._model = lambda: model        
    
    if default is None: 
      try:
        default = form.default
      except AttributeError:
        pass
    self._typedefault = typedefault
    if default is None: default = typedefault
    self._default = default    
    
    self._parent = parent
    self._membername = membername
    assert (parent is None) == (membername is None), (parent, membername)    
    self._toplevel = (parent is None)    
    
    hardform = None #valid value! it means "defined by the parent"
    if typetree is not None and typetree.is_default: 
      hardform = False          
    if hasattr(form, "form") and form.form is not None:
      if form.form not in ("hard", "soft"):
        raise TypeError("Cannot create controller for %s (%s): spyderform 'form' attribute must be soft or hard, not '%s''" % (form.typename, prenam[:-1]), form.form)    
      hardform = True if form.form == "hard" else False
    if self._toplevel: 
      hardform = True
    
    self._hardform = hardform
    self._active = False
    if self._toplevel: self._name = None
    else: self._name = prename
    
    self._views = {}
    self._dynamic_array = False

  def _view(self, viewname=None):
    if viewname is None:
      if len(self._views) == 1: 
        viewname = list(self._views.keys())[0]
      elif "default" in self._views:
        viewname = "default"
      elif len(self._views) == 0:
        raise ValueError("Cannot return view: No views have been defined")        
      else:
        raise ValueError("Cannot return view: there are multiple views, and view name is ambiguous")
    if viewname not in self._views: raise ValueError("Cannot return view: unknown view name \"%s\"" % viewname)  
    return self._views[viewname]
    
  def _is_essential(self):    
    # a child is essential if:
    # - its hardmode is None 
    # - its parent is active and essential
    # - it is a non-optional member of the parent
    # "essential" children behave similar to hardmode 
    # Therefore, is basically a dynamic kind of hardmode
    if self._hardform == True: return True
    if self._hardform == False: return False
    m = self._model()
    if m._is_optional: return False
    if self._parent is None: return True
    return self._parent._active and self._parent._is_essential()
    
  def _listen(self):
    if self._listening: return
    self._listening = True
    #assert self._toplevel or self._parent._model()._type is None
    if self._model()._type is not None:
      self._model()._listen(self._comparator._set)
    else:
      for child in self._children.values():
        child._listen()
  def _unlisten(self):
    if not self._listening: return
    self._listening = False
    #assert self._toplevel or self._parent._model()._type is None
    if self._model()._type is not None:
      self._model()._unlisten(self._comparator._set)
    else:
      for child in self._children.values():
        child._unlisten()
      
        
class controller_leaf(controller_base):
  def __init__(self, 
   form, model, typedefault=None,default=None, 
   prename=None,typetree=None, parent=None, membername=None, model_is_func = False
  ):
    controller_base.__init__(self,
     form,model,typedefault,default,
     prename,typetree,parent,membername,model_is_func
    )  
    self._comparator = comparator(self, True)
    self._formdefault = self._default
    self._value = None
    self._value2 = None
    self.__dif = None
    self.__truedif = None
    self._feedback_cycle = 0
    self._guess_value = None
    self._embedded_resource = False
    if _nonevalue(self._typedefault):
      self._always_active = False
    else:
      self._always_active = True
      if self._typedefault == self._formdefault: self._always_active = False
    self._busy = False    
    
  def _do_set(self, value, extern, forced = False):
    if self._embedded_resource: return
    if value is None: value = self._typedefault
    if not extern: 
      self._model()._set(value)
    if forced: return
    if not len(self._views): return        
    
    if value == self._typedefault: value = None
    v = value    
    if v is None: v = self._formdefault
    if v is None:
      if hasattr(self._form, "options") and self._form.options is not None:
        v = None
      else:
        view = list(self._views.values())[0]
        vval = view.value
        if isinstance(vval, int): v = 0
        elif isinstance(vval, float): v = 0.0
        elif isinstance(vval, str): v = ""
        elif isinstance(vval, bool): v = False
        elif isinstance(vval, tuple) and len(vval) == 2 and \
         (isinstance(vval[0], str) or vval[0] is None) and validvar2(vval[1]):
            v = ""  #file
        else: raise Exception(vval) #cannot guess the null value of the widget
    elif isinstance(v, tuple) and len(v) == 2 and (isinstance(v[0], str) or v[0] is None) and validvar2(v[1]):
      v = v[0]
      if v is None: v = ""
    for viewname, view in self._views.items():  
      view.set(v)
    
  def _set(self, value, extern=False):
    if self._busy:
      if value != self._value_original:
        self._feedback = True
        self._feedbackvalue = value
      return
        
    self._busy = True
    self._feedback = False
    self._value_original = value
    
    if self._model()._is_file or self._model()._is_resource:      
      self._embedded_resource = False
      if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str) and validvar2(value[1]):        
        if value[0] == None: value = None
        elif value[0] == "": value = None        
        elif self._model()._is_resource:
          if value[0] == "<Embedded resource>": 
            self._embedded_resource = True
            value = self._model()._get()
    implicit_none = False
    if not extern: 
      if self._formdefault == value or self._formdefault is None:
        if _get_hardform(self) == False and _nonevalue(value) == True:      
          implicit_none = True
    #special case: options
    if hasattr(self._form, "options") and self._form.options is not None \
     and value is None:
      implicit_none = True
    fvalue = value    
    
    if (fvalue == self._formdefault or implicit_none) and _get_hardform(self) == False: 
      fvalue = None    
    oldvalue = self._value
    tvalue = value
    if tvalue is None: tvalue = self._formdefault
    if tvalue == self._typedefault: tvalue = None        
    if implicit_none: 
      self._value2 = tvalue
      tvalue = None
    self._value = tvalue
        
    if fvalue is None: 
      self._parent._dif_remove(self._membername)    
    else:       
      self._parent._dif_add(self._membername)    
    self.__dif = fvalue
    try:
      self.__dif = self.__dif.dict()
    except:
      pass
    
    if fvalue is None and _get_hardform(self) == False and not self._always_active: 
      if self._is_essential():
        pass
      elif self._active:
        self._parent._child_deactivates(self._membername,extern=extern)
      elif oldvalue != tvalue: 
        self._do_set(tvalue, extern)
      else:
        self._do_set(tvalue, extern,forced=True)
    else:    
      if not self._active:
        self._parent._child_activates(self._membername, extern=extern)
      elif oldvalue != tvalue: 
        self._do_set(tvalue, extern)
      else:
        self._do_set(tvalue, extern, forced=True)
        
    self._busy = False
    if self._feedback:
      self._feedback_cycle += 1
      if self._feedback_cycle == 10:
        raise Exception("Feedback cycle detected in '%s': '%s' feedbacks to '%s'" % \
         (self._name, self._value_original, self._feedbackvalue))
      self._set(self._feedbackvalue,extern=extern)
    self._feedback_cycle = 0
    
  def _activate(self, extern):
    #print("ACTIVATE-LEAF?", self._name)
    # if there is no dif and hardform is False, refuse
    if self.__dif is None and _get_hardform(self) == False and not self._always_active and not self._is_essential(): return
    # else, activate        
    if self._always_active: self._value = self._formdefault
    #if not self._active: print("ACTIVATE-LEAF!", self._name, self._value, self._guess_value)
    forced = False
    if self._typedefault is None and self._value is None:
      forced = True
      #Hmm, we are forced into hardform, but we have to guess what we have to impose onto the model...            
      # these guesses should be temporary: any view-to-model sync should give us the correct values. 
      #  In fact, when this code is called, we are probably in the middle of a sync            
      if self._guess_value is not None:
        #from an earlier view-to-model sync
        self._value = self._guess_value
      elif self._parent._arraycount:
        self._value = None
      elif self._formdefault is not None:
        #the formdefault
        self._value = self._formdefault
      elif self._value2 is not None:
        #a previously imposed "implicit none"
        self._value = self._value2        
      elif hasattr(self._form, "options") and self._form.options is not None:
        #an unselected option
        self._value = None        
      else:  
        # we really have to guess...
        formtype = getattr(self._form, "type", "text")
        if formtype == "number": self._value = 0
        elif formtype in ("checkbox", "switch"): self._value = False
        elif formtype in ("text", "textarea"): self._value = ""
        elif formtype == "file": self._value = None        
        else: self._value = ""
    self._active = True
    self._do_set(self._value,extern,forced=forced)
    
  def _deactivate(self,extern):
    assert self._always_active == False, self._name
    assert self.__dif is None    
    #print("DEACTIVATE-LEAF?", self._name)
    # if hardform is True, refuse
    if _get_hardform(self) == True: return
    #if self._active: print("DEACTIVATE-LEAF!", self._name)
    self._active = False
    self._do_set(None,extern)
    
  def _dif(self):
    return self.__dif    

  def _truedif(self):
    if self._model()._constructionstate == "empty": return None
    return self.__dif    
    
  def _bind_view(self, view, viewname="default"):
    if viewname in self._views: raise ValueError("Cannot bind view: view name \"%s\" already bound" % viewname)
    self._views[viewname] = view
    view.listen(self._set)
  
  def _unbind_view(self, viewname="default"):
    if viewname not in self._views: raise ValueError("Cannot unbind view: unknown view name \"%s\"" % viewname)
    self._views.pop(viewname)
    view.unlisten(self._set)    
    
  def _sync_from_view(self, viewname):   
    if not len(self._views): return
    view = self._views[viewname]
    value = view.value
    if self._guess_value is None:
      self._guess_value = value
    self._set(value)

  def _clear(self):  
    self._active = False  
    self._guess_value = None
    self._value = None
    self._value2 = None
    self.__dif = None
    self._do_set(None,extern=False)
    
class controller(controller_base):
  def __init__(self, 
   form, model, typedefault=None,default=None, 
   prename=None,typetree=None, parent=None, membername=None, model_is_func = False
  ):
    if typetree is None and form is not None: typetree = form._typetree
    controller_base.__init__(self,
     form,model,typedefault,default,
     prename,typetree,parent,membername,model_is_func
    )    
    self._comparator = comparator(self, False)
    prenam = ""
    if prename is not None: prenam = prename + "-"
    self._prenam = prenam

    self._arraycount = form.arraycount
    if self._arraycount > 0 and (not hasattr(form, "form") or form.form is None) :
      raise TypeError("Cannot create controller for %s (%s): spyderform is an array and has no 'form' attribute (soft or hard)'" % (form.typename, prenam[:-1]))    
    
    self._children = {}
    self._difchildren = set()
    self._form = form
    if self._arraycount:
      assert hasattr(form, "length") and form.length is not None
      self._formdefault = {}
      for n in range(form.length):
        self._add_arraychild(n)
    else:  
      membernames = form.get_membernames()
      assert membernames is not None and len(membernames), membernames
      self._formdefault = {}
      for membername, mtypetree in typetree.members:
        def _submodel(attr):
          return getattr(self._model(), attr)
        subform = getattr(form, membername)
        submodel = partial(_submodel, membername)      
        mdefault = None
        if self._default is not None: mdefault = getattr(self._default, membername) 
        child = create_controller(subform, submodel, 
         mtypetree.default, mdefault,
         prename = prenam + membername,
         typetree = mtypetree,
         parent = self, 
         membername = membername,          
         model_is_func = True
        )
        if child is not None:
          self._children[membername] = child
          if child._formdefault is not None:
            self._formdefault[membername] = child._formdefault
    if not len(self._formdefault): 
      self._formdefault = None
    elif self._arraycount:
      self._formdefault = dict2list(self._formdefault, full=True)    
    
    no_typedef = False
    if self._typedefault is None:
      no_typedef = True
    elif self._arraycount and len(self._typedefault) == 0:
      no_typedef = True
      
    if no_typedef:
      self._always_active = False
    else:
      self._always_active = True
      tdic = self._typedefault.dict()
      if tdic == self._formdefault: 
        self._always_active = False
      else:  
        try:
          v = type(self._typedefault)(self._formdefault)
          if v == self._typedefault: self._always_active = False
        except:
          pass
  
  def _add_arraychild(self, n):
    def _ele_model(attr):
      return self._model()[attr]
    try:
      subform = self._form[n]
    except:
      subform = self._form[None]
    submodel = partial(_ele_model, n)
    mdefault = None
    if self._default is not None and len(self._default) > n: mdefault = self._default[n]                
    child = create_controller(subform, submodel, 
      None,mdefault,
      prename = self._prenam + str(n),
      typetree = subform._typetree,
      parent = self, 
      membername = n,          
      model_is_func = True
    )
    self._children[n] = child
    if child._formdefault is not None:
      try:
        self._formdefault[n] = child._formdefault
      except IndexError:
        for nn in range(len(self._formdefault)+1,n):
          self._add_arraychild(nn)
        self._formdefault.append(child._formdefault)  
    if self._listening: child._listen()
    
  def _dif_add(self, childname):
    if len(self._difchildren) == 0:      
      if self._parent is not None:
        self._parent._dif_add(self._membername)
    self._difchildren.add(childname)
    
  def _dif_remove(self, childname):  
    if childname in self._difchildren:
      self._difchildren.remove(childname)
    if len(self._difchildren) == 0:      
      if self._parent is not None:
        self._parent._dif_remove(self._membername)
      
  def _child_activates(self, childname, extern):
    """
    The child indicates to me (the parent) that it wants to activate
    Three possible responses:
    1. I am already activated; in that case, activate the child
    2. The child activates me; in that case, pass on the activation signal to my parent. 
       If I have no parent, I activate myself
    3. I remain inactive even after this; in that case, do nothing (i.e. forbid the child to activate)
    """    
    #print("ACTIVATE-CHILD", self._name, childname)
    if self._active: #1
      self._children[childname]._activate(extern=extern)
      if self._dynamic_array: arraymanager_dynamic.update_visibility(self)      
      return
    if len(self._difchildren) or self._always_active: #2
      if self._toplevel: 
        self._activate(extern=extern)
      else:
        self._parent._child_activates(self._membername,extern=extern)
      if self._dynamic_array: arraymanager_dynamic.update_visibility(self)        
    else: #3
      pass
        
  def _child_deactivates(self, childname, extern):
    """
    The child indicates to me (the parent) that it wants to deactivate
    Three possible responses:
    - I am already deactivated; this can never happen, must be a bug
    - The child deactivates me; in that case, pass on the deactivation signal to the parent
       If I have no parent, I deactivate myself
    - I remain active even after this; in that case, deactivate the child (take care in case of arrays)
    """
    #print("DEACTIVATE-CHILD", self._name, childname)
    assert self._active #1
    if not self._always_active and len(self._difchildren) == 0: #2
      if self._toplevel: 
        self._deactivate(extern=extern)
      else:
        self._parent._child_deactivates(self._membername,extern=extern)
      if self._dynamic_array: arraymanager_dynamic.update_visibility(self)
    else: #3
      self._children[childname]._deactivate(extern=extern)
      if self._arraycount and not _get_hardform(self):
        maxnr = _max_array_child(self)
        for n in range(maxnr, childname):
          self._children[n]._deactivate(extern=extern)
      if self._dynamic_array: arraymanager_dynamic.update_visibility(self)    
          
  def _activate(self, extern):
    #if there is no dif and hardform is False, refuse
    if not len(self._difchildren) and _get_hardform(self) == False and not self._is_essential(): return
    #else, activate
    self._active = True
    if self._arraycount and not _get_hardform(self):
      maxnr = _max_array_child(self)
      for n in range(maxnr):
        child = self._children[n]
        child._activate(extern=extern)
    else:  
      for childname in sorted(self._children):
        child = self._children[childname]
        child._activate(extern=extern)
    if self._dynamic_array:
      arraymanager_dynamic.update_visibility(self)
      
  def _deactivate(self, extern):
    assert len(self._difchildren) == 0
    assert self._always_active == False, self._name
    if self._toplevel:
      self._model()._set(None)    
    #if hardform is True, refuse
    if _get_hardform(self) == True: return
    #print("DEACTIVATE", self._name)
    self._active = False
    for childname in sorted(self._children):
      child = self._children[childname]
      child._deactivate(extern=extern)
    if self._dynamic_array:
      arraymanager_dynamic.update_visibility(self)

  def __len__(self):
    assert self._arraycount > 0      
    return max(self._children.keys())+1
  def __getitem__(self, index):
    assert isinstance(index, int)
    assert self._arraycount > 0
    assert index >= 0 #no negative indices
    if index not in self._children:
      for n in range(max(self._children)+1,index+1):
        self._add_arraychild(n)        
    return self._children[index]
  def __setitem__(self, index, value):
    assert isinstance(index, int)
    assert self._arraycount > 0
    assert index >= 0 #no negative indices
    self._children[index]._set(value)
  def __getattr__(self, attr):  
    if attr.startswith("_"): raise AttributeError(attr)
    try:
      return self._children[attr]
    except KeyError:
      raise AttributeError(attr)
  def __setattr__(self, attr, value):
    if attr.startswith("_"): 
      self.__dict__[attr] = value
    else:  
      self._children[attr]._set(value)

  def _dif(self): 
    if not len(self._difchildren): return None
    ret = {}
    for child in self._difchildren:
      ret[child] = self._children[child]._dif()
    if self._arraycount:      
      ret = dict2list(ret,full=True)
    return ret  
      
  def _truedif(self): 
    if not len(self._difchildren): return None
    m = self._model()
    if m._constructionstate == "empty": return None
    if not self._toplevel:
      p = self._parent
      mpd = p._model()._model__ddif    
      if mpd is None or m._membername not in mpd: return None
        
    ret = {}
    for child in self._difchildren:
      d = self._children[child]._truedif()
      if d is not None: ret[child] = d
    if not len(ret): return None
    if self._arraycount:      
      ret = dict2list(ret,full=True)
    return ret  
    
  def _bind_view(self, view, viewname="default"):
    if viewname in self._views: raise ValueError("Cannot bind view: view name \"%s\" already bound" % viewname)
    self._views[viewname] = view
    for childname in sorted(self._children):
      child = self._children[childname]
      if self._arraycount:
        childview = view[childname]
      else:
        try:
          childview = getattr(view,childname)
        except AttributeError:
          continue #no editable widgets for this child...
      child._bind_view(childview, viewname)
    if len(self._views) == 1:
      if self._parent is None:
        arraymanager_dynamic.configure(self) #sets "_dynamic_array" attribute on dynamic-array subforms
        arraymanager_dynamic.update_all_visibility(self)

  
  def _unbind_view(self, viewname="default"):
    if viewname not in self._views: raise ValueError("Cannot unbind view: unknown view name \"%s\"" % viewname)
    self._views.pop(viewname)
    for childname in sorted(self._children):
      child = self._children[childname]
      child._unbind_view(viewname)
      
  def _sync_to_view(self):
    assert self._toplevel
    value = self._model()._truevalue
    self._comparator._set(value)
 
  def _sync_from_view(self, viewname=None):
    m = self._model()
    m._disable_passdown()
    m._start_syncing()
    if viewname is None:
      if len(self._views) == 1: 
        viewname = list(self._views.keys())[0]
      elif "default" in self._views:
        viewname = "default"
      else:
        raise ValueError("Cannot sync from view: there are multiple views, and view name is ambiguous")
    if viewname not in self._views: raise ValueError("Cannot sync from view: unknown view name \"%s\"" % viewname)
    
    for childname in sorted(self._children):
      child = self._children[childname]
      child._sync_from_view(viewname)
    m._enable_passdown()    
    m._stop_syncing()
    
  def _clear(self):  
    self._active = False    
    m = self._model()
    m._disable_passdown()    
    for child in self._children.values():
      child._clear()
    m._enable_passdown()
    
  def _insert(self, index):
    assert self._arraycount
    self._model()._insert(index)
    difchildren = set()
    actions = []
    for n in range(index-1):
      if n in self._difchildren:
        difchildren.add(n)
    for n in reversed(range(index, self._form.length-1)):
      f2 = self._children[n+1]
      f1 = self._children[n]
      if n in self._difchildren: difchildren.add(n+1)
      actions += copy_controller(f1, f2)      
    self._children[index]._clear()  
    for a,v in actions: a(v)
    self._difchildren = difchildren  

  def _delete(self, index):
    assert self._arraycount
    assert index >= 0 and index < self._form.length, index
    ok = self._model()._delete(index)
    if ok:
      l = len(self._model()._get())
      self._children[l]._clear()
    else:
      difchildren = set()
      actions = []
      l = self._form.length-1
      for n in range(index-1):
        if n in self._difchildren:
          difchildren.add(n)
      for n in range(index, l):
        f2 = self._children[n+1]
        f1 = self._children[n]
        if (n+1) in self._difchildren: difchildren.add(n)
        actions += copy_controller(f2, f1)
      if l in self._difchildren: 
        self._children[l]._clear()
      for a,v in actions: a(v)  
      self._difchildren = difchildren  
    
    
