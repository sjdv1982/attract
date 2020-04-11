from .. import Handler
from ..Handler import Method, PassDown, NullHandler
from ..hamlet.tag import SimpleTag


def get_tip(node):
  ch = node.getChildren()
  if len(ch) == 0: return node
  elif len(ch) > 1: 
    if node.space.htmltip is not None:
      return get_tip(node.space.htmltip)
    else:
      return None
  only_child = node[ch[0]]
  return get_tip(only_child)
    
def morph_markup(node):
  space = node.space
  assert node.space.htmlnode is not None
  htmlnode = node.space.htmlnode
  tip_h = node.space.htmltip
  if tip_h is None:
    tip_h = get_tip(htmlnode)  
  if tip_h is None and node.getChildren():
    raise TypeError("html-concrete-node %s\n:  cannot make markup, both node and htmlnode have multiple children" % str(space.formpath))
  if not node.getChildren(): 
    #replace node with htmlnode
    targetspace = htmlnode.space
    for var in vars(space):
      if var == "htmlnode": continue
      setattr(targetspace, var, getattr(space, var))
    assert node.root is not None  
    node.parent.replaceChild(node.name, htmlnode)  
  elif tip_h is not htmlnode:  
    #morph node into tip
    Handler.passdown("make-markup", node)
    sourcespace = tip_h.space
    for var in vars(sourcespace):
      setattr(space, var, getattr(sourcespace, var))
    assert node.root is not None, (node.space.formpath, node.nodetype )
    node.morph(tip_h.nodetype)    
    #attach htmlnode in the place of node
    node.parent.detachAndReplaceChild(node.name, htmlnode)
    #attach node in the place of tip    
    assert tip_h.parent is not None
    tip_h.parent.replaceChild(tip_h.name, node)    
  else: 
    #morph node into htmlnode
    Handler.passdown("make-markup", node)
    sourcespace = htmlnode.space
    for var in vars(sourcespace):
      setattr(space, var, getattr(sourcespace, var))
    assert node.root is not None, (node.space.formpath, node.nodetype )
    node.morph(htmlnode.nodetype)
  return Handler.CLAIM

def morph_formcontainer(node):
  if node.parent is not node.root: #only if top element
    return Handler.NO_MATCH 
  Handler.passdown("html-make-abstract", node)  
  node.morph("html-abstract-formcontainer")
  return Handler.CLAIM
      
def init(root):
  root.addHandler( PassDown("html-concrete-node", "html-make-concrete") )
  root.addHandler( Method("html-concrete-node", "make-markup", morph_markup) )
  
  root.addHandler( PassDown("markup-mixedtag", "html-make-concrete") )
  root.addHandler( PassDown("markup-nulltag", "html-make-concrete") )
  root.addHandler( NullHandler("markup-stringtag", "html-make-concrete") )

from . import abstract