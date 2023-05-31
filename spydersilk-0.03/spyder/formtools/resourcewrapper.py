embedded_resource_filename = "<Embedded resource>"

try:
  unicode
except NameError:
  unicode = str

class ResourceWrapper(object):
  def __init__(self, getprop, setprop, typename, formtype):
    import spyder, Spyder
    assert isinstance(typename, str)
    assert typename.startswith("Resource")
    assert formtype in ("file", "text"), formtype #cannot be "data"
    self.resourcetype = getattr(Spyder, typename)
    self.type = getattr(Spyder, typename[len("Resource"):])
    self.formtype = formtype
    self._getprop = getprop
    self._setprop = setprop
    if formtype == "file":
      self.getprop = self._getprop_file
      self.setprop = self._setprop_file
    elif formtype == "text":  
      self.getprop = self._getprop_text
      self.setprop = self._setprop_text
  def _getprop_file(self):
    v = self._getprop()
    if v is None: return None
    if v == embedded_resource_filename: return v
    assert isinstance(v, tuple) and len(v) == 2, v #widget must provide a (filename, filetype) tuple
    return v
  def _setprop_file(self, v):  
    import spyder, Spyder    
    if isinstance(v, str) or isinstance(v, bytes) or isinstance(v, unicode):
      filename = Spyder.String(v)
    elif isinstance(v, Spyder.File):
      filename = v.name
    elif isinstance(v, Spyder.Resource):
      filename = v.filename      
      if filename is None: filename = embedded_resource_filename 
    elif isinstance(v, self.type):
      filename = embedded_resource_filename
    elif v is None:
      filename = None
    else: 
      raise TypeError(type(v))
    self._setprop(filename)
  def _getprop_text(self):
    v = self._getprop()
    if v is None: return None
    val = self.type(v)
    value = self.resourcetype(data=val)
    if value._data is None: return None
    return value
  def _setprop_text(self, v):
    if v is None: 
      self._setprop(None)  
    else:
      value = self.resourcetype(data=v)
      self._setprop(str(value.data()))  
    