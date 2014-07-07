from spyder.stalkwalk import Handler
from spyder.stalkwalk.Handler import Method, CMethod, FallBack, NullHandler, PassDown
from spyder.stalkwalk.hamlet.tag import SimpleTag, VoidTag, ContainerTag, NullTag
from .concrete import *
  
def morph_inputleaf(node):    
  #if node.space.complexity > 0: return Handler.NO_MATCH #can never happen...
  space = node.space
  if space.formtype in ("none", None): return Handler.NO_MATCH
  node.morph("html-abstract-inputleaf")
  return Handler.CLAIM


def morph_blockcategory(node):  
  g0 = node.space.group
  if g0 is None or g0.type != "elegroup": return Handler.NO_MATCH  
  g = g0.elegroup  
  if g is None: return Handler.NO_MATCH  
  if g.type != "category": return Handler.NO_MATCH 
  if node.space.complexity > 1: return Handler.NO_MATCH 
  assert node.space.wrapping == 1, (node.name, node.space.formpath, node.space.wrapping)  
  
  node.space.blockname = g.categoryname
  if hasattr(g, "controltitle"):
    node.space.controltitle = g.controltitle
  else:
    node.space.controltitle = g.title
  
  x = "temporary-category-childname"
  node.addChild("html-abstract-category", x)
  c = node.detachChild(x)
  c.space.blockid = g.id
  c.space.categoryname = g.categoryname
  c.space.page = g.page
  c.space.icon = g.icon
  c.space.title = g.title
  if hasattr(g, "title2"):
    c.space.title2 = g.title2
  c.space.description = g.description
  if hasattr(g, "html_description"):
    c.space.html_description = g.html_description
      
  for var in vars(node.space):
    setattr(c.space, var, getattr(node.space, var))
  node.parent.detachAndReplaceChild(node.name, c)
  c.attachChild(node, "block")
  node.morph("html-abstract-level1-block")
  Handler.passdown("html-make-abstract", node)  
  return Handler.CLAIM

def morph_category(node):
  g0 = node.space.group
  if g0 is None or g0.type != "elegroup": return Handler.NO_MATCH  
  g = g0.elegroup  
  if g is None: return Handler.NO_MATCH  
  if g.type != "category": return Handler.NO_MATCH 
  if node.space.complexity == 1: return Handler.NO_MATCH 
  assert node.space.wrapping == 1, (node.name, node.space.formpath, node.space.wrapping)  
  
  node.space.categoryid = g.id
  node.space.categoryname = g.categoryname
  node.space.page = g.page
  node.space.icon = g.icon
  node.space.title = g.title
  if hasattr(g, "title2"):
    node.space.title2 = g.title2  
  node.space.description = g.description
  if hasattr(g, "html_description"):
    node.space.html_description = g.html_description    
    
  node.morph("html-abstract-category")
  Handler.passdown("html-make-abstract", node)
  return Handler.CLAIM

def morph_block(node):
  g0 = node.space.group
  if g0 is None or g0.type != "elegroup": return Handler.NO_MATCH  
  g = g0.elegroup  
  if g is None: return Handler.NO_MATCH  
  if g.type != "block": return Handler.NO_MATCH    
  
  w = node.space.wrapping
  if w == 2:
    node.morph("html-abstract-level1-block")
  elif w == 3 or (w == 4 and node.parent.nodetype == "html-abstract-level1-cloneblock"):
    node.morph("html-abstract-level2-block")    
  elif w == 4 or (w == 5 and node.parent.parent.nodetype == "html-abstract-level1-cloneblock"):
    node.morph("html-abstract-level3-block")    
  g = node.space.group.elegroup  
  
  if hasattr(g, "blockname"):
    node.space.blockname = g.blockname
  else:
    node.space.blockname = g.id
  if hasattr(g, "title"):
    node.space.title = g.title 
  if hasattr(g, "title2"):
    node.space.title2 = g.title2
  if hasattr(g, "controltitle"):
    node.space.controltitle = g.controltitle
  elif hasattr(g, "title"):
    node.space.controltitle = g.title
  
  Handler.passdown("html-make-abstract", node)  
  return Handler.CLAIM
  
def morph_fallback(node):
  if node.parent and node.parent is node.root: #top form element
    return Handler.NO_MATCH   
  if node.space.wrapping == 1: return Handler.NO_MATCH  #can't make categories
  if node.space.wrapping == 2: return Handler.NO_MATCH  #can't make blocks
  
  w = node.space.wrapping
  if w == 3 or (w == 4 and node.parent.nodetype == "html-abstract-level1-cloneblock"):
    assert node.space.group is None or node.space.group.elegroup is None #else: html-abstract-level2-block
    node.morph("html-abstract-level2")    
  elif w == 4 or (w == 5 and node.parent.parent.nodetype == "html-abstract-level1-cloneblock"):
    node.morph("html-abstract-level3")    
  else:  
    raise Exception(node.name, node.space.formpath, node.space.wrapping) #too much wrapping...
    
  Handler.passdown("html-make-abstract", node)    
  return Handler.CLAIM
  
def morph_cloneblock(node):
  if node.parent is None: return Handler.NO_MATCH
  if node.parent.nodetype != "html-abstract-clonecontainer": return Handler.NO_MATCH
  node.morph("html-abstract-level1-cloneblock")
  Handler.passdown("html-make-abstract", node)  
  return Handler.CLAIM
  
def morph_clonecontainer(node):
  from .concrete import _extract_fields
  
  space = node.space
  if space.form is None or space.form.arraycount == 0: return Handler.NO_MATCH  

  fields = ("clonebutton","clonelength","controltitle", "title", "blockname")
  p = _extract_fields(space, space.form, fields)  
  if p.clonelength is None: return Handler.NO_MATCH
  assert node.parent.nodetype == "html-abstract-category", (node.name, node.formpath) #only categories can clone
    
  clone = Space.Space(top=False)
  clone.button = p.clonebutton
  assert p.clonelength is not None, node.name
  clone.length = p.clonelength
  if p.blockname is None: p.blockname = node.name
  clone.name = p.blockname
  
  node.morph("html-abstract-clonecontainer")
  Handler.passdown("html-make-abstract", node)  
  if node.parent.space.clones is None: node.parent.space.clones = []
  node.parent.space.clones.append(clone)
  children = node.getChildren()
  for n in children:
    child = node[children[n]]
    child.space.clone = clone
    child.space.blockindex = n
    if p.controltitle is None: p.controltitle = p.title
    if p.controltitle is None: p.controltitle = p.name
    child.space.controltitle = p.controltitle + " " + str(n+1) 
    child.space.blockname = p.blockname
    
  return Handler.CLAIM

def morph_formcontainer2(node):
  if node.parent is not node.root: #only if top element
    return Handler.NO_MATCH   
  from spyder.stalkwalk.html import morph_formcontainer
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

def morph_nothing(node):    
  space = node.space
  if space.formtype in ("none", None):
    node.parent.removeChild(node.name)
    return Handler.CLAIM
  else:  
    return Handler.NO_MATCH 

################################  
  
def insert_level2(node):    
  Handler.passdown("html-insert-level2", node)
  
  groups = []
  curr_group = None
  for childname in node.getChildren():
    t = node[childname].nodetype
    if t in ("html-abstract-level3", "html-abstract-level3-block", "html-abstract-inputleaf"):
      if curr_group is None: curr_group = []
      curr_group.append(childname)
    else:
      if curr_group is not None: groups.append(curr_group)
      curr_group = None
  if curr_group is not None: groups.append(curr_group)      
  
  if groups:      
    for group in groups:
      position = list(node.getChildren()).index(group[0])
      g = node.insertChild("html-abstract-level2", "group-" + group[0], position)
      g.space.complexity = node.space.complexity
      for m in group:
        child = node.detachChild(m)
        g.attachChild(child, m)        
  return Handler.CLAIM

def promote_level3_switch(node):
  #level3 switch blocks are promoted to level2
  has_switch = False
  if node.space.group and node.space.group.elegroup:
    if hasattr(node.space.group.elegroup, "has_switch"):
      has_switch = node.space.group.elegroup.has_switch
  if not has_switch: return
  
  node.morph("html-abstract-level2-block")
  granpa = node.parent.parent
  pname = node.parent.name
  pos = granpa.getChildren().index(pname)
  
  node.parent.detachChild(node.name)
  granpa.attachChild(node,pname + "-" + node.space.group.name,pos+1)
  
def _init_concrete(root):
  root.addHandler( Method("html-abstract-level2", "html-make-concrete", morph_level2) )
  root.addHandler( Method("html-abstract-level3", "html-make-concrete", morph_level3) )
  root.addHandler( Method("html-abstract-category", "html-make-concrete", morph_categorywrapper) )
  root.addHandler( Method("html-abstract-clonecontainer", "html-make-concrete", morph_clonewrapper) )
  root.addHandler( Method("html-abstract-level1-cloneblock", "html-make-concrete", morph_cloneblockwrapper) )
  root.addHandler( Method("html-abstract-level1-block", "html-make-concrete", morph_level1_blockwrapper) )
  root.addHandler( Method("html-abstract-level2-block", "html-make-concrete", morph_level2) )
  root.addHandler( Method("html-abstract-level3-block", "html-make-concrete", morph_level3) )
  root.addHandler( Method("html-abstract-formcontainer", "html-make-concrete", morph_formwrapper) )
  
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_text) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_number) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_textarea) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_selecttag) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_file) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_checkbox) )
  root.addHandler( CMethod("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_switch) )
  root.addHandler( FallBack("html-abstract-inputleaf", "html-make-concrete", morph_inputtag_unknown) )
  
  root.addHandler( Method("html-concrete-fieldcontainer", "make-markup", morph_markup_fieldcontainer) )
  root.addHandler( NullHandler("html-concrete-fieldcontainer", "html-make-concrete") )
  
  root.addHandler( NullHandler("html-abstract-level2", "html-insert-level2") )
  root.addHandler( NullHandler("html-abstract-level3", "html-insert-level2") )
  root.addHandler( PassDown("html-abstract-category", "html-insert-level2") )
  root.addHandler( PassDown("html-abstract-clonecontainer", "html-insert-level2") )
  root.addHandler( Method("html-abstract-level1-cloneblock", "html-insert-level2", insert_level2) )
  root.addHandler( Method("html-abstract-level1-block", "html-insert-level2", insert_level2) )
  root.addHandler( NullHandler("html-abstract-level2-block", "html-insert-level2") )
  root.addHandler( NullHandler("html-abstract-level3-block", "html-insert-level2") )
  root.addHandler( PassDown("html-abstract-formcontainer", "html-insert-level2") )
  root.addHandler( NullHandler("html-abstract-inputleaf", "html-insert-level2") )

  root.addHandler( PassDown("html-abstract-level2", "html-promote-level3-switch") )
  root.addHandler( NullHandler("html-abstract-level3", "html-promote-level3-switch") )
  root.addHandler( PassDown("html-abstract-category", "html-promote-level3-switch") )
  root.addHandler( PassDown("html-abstract-clonecontainer", "html-promote-level3-switch") )
  root.addHandler( PassDown("html-abstract-level1-cloneblock", "html-promote-level3-switch") )
  root.addHandler( PassDown("html-abstract-level1-block", "html-promote-level3-switch") )
  root.addHandler( PassDown("html-abstract-level2-block", "html-promote-level3-switch") )
  root.addHandler( Method("html-abstract-level3-block", "html-promote-level3-switch", promote_level3_switch) )
  root.addHandler( PassDown("html-abstract-formcontainer", "html-promote-level3-switch") )
  root.addHandler( NullHandler("html-abstract-inputleaf", "html-promote-level3-switch") )
  
  from spyder.stalkwalk.html import init
  init(root)
  
def init_spyderform(root):
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_fallback) )  
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_formcontainer2) )
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_clonecontainer) )
  root.addHandler( CMethod("spyderform", "html-make-abstract", morph_cloneblock) )

  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_fallback) )  
  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_category) )
  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_blockcategory) )
  root.addHandler( CMethod("spydergroup", "html-make-abstract", morph_block) )
  
  root.addHandler( CMethod("spyderwidget", "html-make-abstract", morph_inputleaf) )
  root.addHandler( CMethod("spyderwidget", "html-make-abstract", morph_nothing) )
  _init_concrete(root)
  