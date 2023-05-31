"""
TODO: log access to Space
"""
class Space(object):
  def __init__(self, d={}, top=True):
    for key,value in d.items():
      setattr(self, key, value)
    if top:
      self.params = Space(top=False)
  def __setattr__(self, attr, value):
    if value is None and attr in self.__dict__:
      self.__dict__.pop(attr)
    elif value is not None:
      self.__dict__[attr] = value
  def __getattr__(self, attr):    
    if attr.startswith("_"): raise AttributeError(attr)
    return None
  def clear(self):
    params = self.params
    self.__dict__.clear()
    self.params = params
  def __str__(self):
    return str(self.__dict__)
    
def logger(*args):
  #to be overridden...
  pass
  
def _logger(*args):
  return logger(*args)