import spyder, Spyder
class comparator(object):
  def __init__(self, controller, leaf):
    self._value = None
    self._busy = False
    self._defined_value = False   
    self._strvalue = None
    
    assert leaf in (True, False), leaf
    self._controller = controller
    self._leaf = leaf
  def _set(self, value):
    if self._defined_value and value == self._value: return
    self._value = value
    self._defined_value = True
    
    con = self._controller    
    if self._leaf:
      con._set(value, extern=True)
      return

    assert not self._busy, self._controller._name
    strvalue = str(value)
    if strvalue == self._strvalue: return
    self._strvalue = strvalue
    self._busy = True      
    
    model = con._model()
    if isinstance(value, Spyder.Resource):
      value = value.data()
    if value is None:
      for childname in sorted(con._children):
        child = con._children[childname]
        child._comparator._set(None)
    elif con._arraycount:
      for childname in sorted(con._children):        
        try:
          subvalue = value[childname]
        except IndexError:
          subvalue = None
        child = con._children[childname]
        child._comparator._set(subvalue)
    elif model._type is None:
      for childname in sorted(con._children):        
        try:
          subvalue = value[childname]
        except KeyError:
          subvalue = None
        child = con._children[childname]
        child._comparator._set(subvalue)        
    else:
      for childname in sorted(con._children):        
        subvalue = None
        try:
          subvalue = getattr(value, childname)
        except AttributeError:
          subvalue = None
        child = con._children[childname]
        child._comparator._set(subvalue)
    self._busy = False    
      