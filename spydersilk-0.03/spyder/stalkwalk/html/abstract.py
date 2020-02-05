from .. import Handler
from ..Handler import Method, CMethod, FallBack
from ..hamlet.tag import SimpleTag, VoidTag, ContainerTag, NullTag
from .concrete import *
  
def morph_inputleaf(node):    
  #if node.space.complexity > 0: return Handler.NO_MATCH #can never happen...
  space = node.space
  if space.formtype in ("none", None): return Handler.NO_MATCH
  node.morph("html-abstract-inputleaf")
  return Handler.CLAIM

def morph_h1(node):
  if node.space.wrapping == 1: return Handler.NO_MATCH
  if node.space.complexity > 2: return Handler.NO_MATCH
  if node.space.complexity == 2 and \
   (node.space.form is None or node.space.form.arraycount == 0): return Handler.NO_MATCH
  if node.parent and node.parent is node.root: #top form element
    return Handler.NO_MATCH 
  Handler.passdown("html-make-abstract", node)  
  node.morph("html-abstract-h1")
  return Handler.CLAIM
  
def morph_foldcontainer(node):
  if node.space.complexity == 1 and node.space.wrapping > 1: return Handler.NO_MATCH
  if node.space.complexity == 2 and \
   (node.space.form is not None and node.space.form.arraycount > 0): return Handler.NO_MATCH
  if node.parent and node.parent is node.root: #top form element
    return Handler.NO_MATCH 
  Handler.passdown("html-make-abstract", node)  
  node.morph("html-abstract-foldcontainer")
  return Handler.CLAIM

def morph_formcontainer2(node):
  if node.parent is not node.root: #only if top element
    return Handler.NO_MATCH   
  from . import morph_formcontainer
  morph_formcontainer(node)
  tablenode = SimpleTag("table")
  node.space.htmlnode = tablenode
  tail = ContainerTag(
    "tr", 
    (
      ComplexTag(
       "td", 
       attributes = [
        ("colspan", "2",),
        ("align", "center",),
        ("style", "width:710px"),
       ], 
       children = [
         VoidTag(
          "input",
          [ ("type","submit",{"space_after":1}) ],
          closetag = "slash"
         )
       ],
       inline = True
      ),      
    ),
    lines_before = 2
  )
  tip = NullTag()
  tablenode.attachChild(tip, "tip")
  tablenode.attachChild(tail, "tail")
  tablenode.space.htmltip = tip
  return Handler.CLAIM

def morph_divcontainer(node):
  #if node.space.complexity == 0: return Handler.NO_MATCH #can never happen...
  if node.parent and node.parent is node.root: #top form element
    return Handler.NO_MATCH 
  Handler.passdown("html-make-abstract", node)  
  node.morph("html-abstract-divcontainer")
  return Handler.CLAIM

def morph_nothing(node):    
  space = node.space
  if space.formtype in ("none", None):
    node.parent.removeChild(node.name)
    return Handler.CLAIM
  else:  
    return Handler.NO_MATCH 
  
def _init_concrete(root):
  root.addHandler( Method("html-abstract-divcontainer", "html-make-concrete", morph_divwrapper) )
  root.addHandler( Method("html-abstract-h1", "html-make-concrete", morph_h1wrapper) )
  root.addHandler( Method("html-abstract-foldcontainer", "html-make-concrete", morph_foldwrapper) )
  root.addHandler( Method("html-abstract-formcontainer", "html-make-concrete", morph_formwrapper) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_text) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_text_spin) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_textarea) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_selecttag) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_file) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_checkbox) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_hidden) )
  root.addHandler( FallBack("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_unknown) )
  from . import init
  init(root)
  
def init_spyderform(root):
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_divcontainer) )  
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_formcontainer2) )
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_foldcontainer) )
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_h1) )  

  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_divcontainer) )
  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_foldcontainer) )
  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_h1) )
  
  root.addHandler( CMethod("spyderwidget", "html-make-abstract", morph_inputleaf) )
  root.addHandler( CMethod("spyderwidget", "html-make-abstract", morph_nothing) )
  _init_concrete(root)
  