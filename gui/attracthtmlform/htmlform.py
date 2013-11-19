from spyder.stalkwalk import spyderform
from spyder.stalkwalk import html
from . import abstract
import spyder
from Spyder import Object
class HtmlformRoot(spyderform.SpyderformRoot): pass

def _htmlform(
 root,obj,form,header,footer, 
 indent, header_indentation, 
 cgi, hidden
):     
  from spyder.stalkwalk.testing import print_grouping, print_tree
  from spyder.stalkwalk import markup
  #print("ORIGINAL SPYDERFORM")
  #print_grouping(root)  
  #print
  #print_tree(root)
  
  markup.init(root)  
  #print("HTML-MAKE-ABSTRACT")
  root["top"].space.cgi = cgi
  root["top"].space.hidden = hidden
  root["top"].sendMessage("html-make-abstract")
  root["top"].sendMessage("html-insert-level2")
  root["top"].sendMessage("html-promote-level3-switch")
  #print_tree(root)
  #import sys; sys.exit()
  
  #print("HTML-MAKE-CONCRETE")
  root["top"].sendMessage("html-make-concrete")
  #print_tree(root)
  #import sys; sys.exit()
  
  #print("HTML-MAKE-MARKUP")
  root["top"].sendMessage("make-markup")  
  #print_tree(root)
  
  #print("HTML-COMPLETE-MARKUP")
  root["top"].sendMessage("complete-markup")
  #print_tree(root)
  
  root.space.indent = indent
  root.space.printer.clear()
  root.space.printer.value = header
  root.space.printer.indentation = header_indentation
  #print("HTML-PRINT-MARKUP")
  root["top"].sendMessage("print-markup")
  htmlcode = root.space.printer.value
  htmlcode += footer
  return htmlcode


def htmlform(
 obj=None, form=None,
 header = "", footer = "", 
 indent = 2, header_indentation = 0,
 cgi = None, hidden = None
): 
  """
  obj can be: A Spyder class, a Spyder object, or None
    if obj is None, then form must not be None
  form is a spyderform object
    if form is None, it is extracted from the Spyder obj
  """
  assert cgi is not None
  assert isinstance(hidden, dict) or hidden is None, hidden
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
  abstract.init_spyderform(root)
  return _htmlform(
   root,obj,form,header,footer, 
   indent, header_indentation, 
   cgi, hidden
  ) 
  