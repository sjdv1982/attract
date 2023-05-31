import sys, os
class moduleSpyder_wrapper(object):
  __path__ = os.path.split(__file__)[0]
  def __getattr__(self, attr):
    import spyder, spyder.spydercompile
    if attr == "__all__":
      dic = spyder.load(attr)
      for d in dic: globals()[d] = dic[d]
      return list(dic.keys())
    elif spyder.spydercompile.validvar2(attr) or attr == "Object":
      return spyder.load(attr)
    try:
      ret = globals()[attr] #non-capital: no wrapping at all
    except KeyError:
      raise AttributeError(attr)

sys.modules["Spyder"] = moduleSpyder_wrapper()
