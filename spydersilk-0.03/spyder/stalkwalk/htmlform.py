from . import spyderform
from . import markup
from . import html
class HtmlformRoot(spyderform.SpyderformRoot): pass

import spyder
from Spyder import Object

def _htmlform(
 root,obj,form,header,footer, 
 indent, header_indentation, 
 cgi, hidden, resourcefilename, newtab
):     
  markup.init(root)  
  if resourcefilename is not None:
    assert form.resourcefilevar is not None
    hidden[form.resourcefilevar] = resourcefilename
  root.space.arraymarker = getattr(form, "arraymarker", None)
  root["top"].space.cgi = cgi  
  root["top"].space.hidden = hidden
  root["top"].space.newtab = newtab
  root["top"].sendMessage("html-make-abstract")
  root["top"].sendMessage("html-make-concrete")
  root["top"].sendMessage("make-markup")
  root["top"].sendMessage("complete-markup")
  
  root.space.indent = indent
  root.space.printer.clear()
  root.space.printer.value = header
  root.space.printer.indentation = header_indentation
  root["top"].sendMessage("print-markup")
  htmlcode = root.space.printer.value
  htmlcode += footer
  return htmlcode

def htmlform(
 obj=None, form=None,
 header = "", footer = "", 
 indent = 2, header_indentation = 0,
 cgi = None, hidden = None, resourcefilename = None, 
 newtab = False
): 
  """
  obj can be: A Spyder class, a Spyder object, or None
    if obj is None, then form must not be None
  form is a spyderform object
    if form is None, it is extracted from the Spyder obj
  """
  assert cgi is not None
  assert isinstance(hidden, dict) or hidden is None, hidden
  if hidden is None: hidden = {}
  else: hidden = hidden.copy()  
  if form is None:
    if isinstance(obj, Object):
      form = obj._form()
    else:
      ok = False
      try: 
        if issubclass(obj, Object):
          form = obj._form() 
          obj = None
          ok = True
      except TypeError:
        pass
      if not ok: raise TypeError(obj)
  else:
    try:
      if issubclass(obj, Object): 
        obj = None
    except TypeError:
      pass    
  assert isinstance(form, spyder.core.spyderform)
  
  root = HtmlformRoot()
  spyderform.build(root, form, obj)
  html.abstract.init_spyderform(root)
  return _htmlform(
   root,obj,form,header,footer, 
   indent, header_indentation, 
   cgi, hidden, resourcefilename, newtab
  ) 
  