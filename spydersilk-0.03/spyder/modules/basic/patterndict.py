"""
A patterndict is a dict that can generate unused keys based on a pattern
"""
class patterndict(dict):
  def __init__(self, defaultpattern, separator="-"):
    self.defaultpattern = defaultpattern
    self.separator = separator
  def new_key(self, pattern = None):
    if pattern == None: pattern=self.defaultpattern
    if pattern not in self: return pattern
    for n in xrange(1,100000):
      key = pattern + self.separator + str(n)
      if key not in self: return key
  
