def _extract_fields(space, form, variables):
  ret = Space.Space(top=False)
  for v in variables:
    if getattr(space, v) is not None:
      setattr(ret, v, getattr(space, v))
    elif form is not None and getattr(form, v) is not None:
      setattr(ret, v, getattr(form, v))
  return ret    

from .concrete_leaf import *
from spyder.stalkwalk import Handler, Space
from spyder.stalkwalk.hamlet.tag import SimpleTag, VoidTag, ComplexTag, ComplexVoidTag, TextWrapTag, StringTag, NullTag

def _add_head(node, tag):
  space = node.space
  formname = space.formname
  if space.formnode is not None:
    formname = space.group.name
  assert formname is not None, (node.name, space.formpath, node.parent.name, node.parent.space.formpath)
  fields = ("tooltip","tooltip_doc","style")
  p = _extract_fields(space, space.form, fields)
  if p.style is None:
    if space.formtype == "switch": 
      p.style = "h3"
    else:
      p.style = "p"
  divtag = SimpleTag ("div",[("class","title")], inline = True)
  style = TextWrapTag(p.style, str(formname), inline = True)
  divtag.attachChild(style, "child")
  tag.attachChild(divtag, "tag-formname")
  if p.tooltip is not None:
    if p.tooltip_doc is not None:
      tooltiptext = TextWrapTag (
       "a",
       p.tooltip,
       [("href",p.tooltip_doc),("target","_blank")],
       inline = True
      )
    else:
      tooltiptext = StringTag(p.tooltip)
    divtag = SimpleTag("div",[("class","tooltip")], inline = True)
    divtag2 = SimpleTag("div",[("class","tooltip-body")], inline = True)
    divtag2.attachChild(tooltiptext, "child")
    divtag.attachChild(divtag2, "child")    
    tag.attachChild(divtag, "tag-tooltip")
  tip = NullTag() 
  tag.attachChild(tip, "tag-tip")
  tag.space.htmltip = tip 
  
def _add_title(node, tag):
  space = node.space
  fields = ("title","title2")
  p = _extract_fields(space, space.form, fields)
  if p.title is None and space.form is not None and node.space.complexity is not None:
    p.title = str(space.formname)
  if p.title is not None: 
    if p.title2 is not None:
      fixtag = SimpleTag("div", [("class","group-inline clearfix")] , 
       comment = " <!--END: group-inline -->", lines_after = 1, lines_inside_after = 1)
       
      divtag = SimpleTag ("div",[("class","field-container")], lines_before = 1)
      divtag2 = SimpleTag ("div",[("class","title")], inline = True)
      h3 = TextWrapTag("h3", p.title, inline = True)
      divtag2.attachChild(h3, "child")
      divtag.attachChild(divtag2, "child")
      fixtag.attachChild(divtag, "tag-title1")
      
      divtag = SimpleTag ("div",[("class","field-container")], lines_before = 1)
      divtag2 = SimpleTag ("div",[("class","title")], inline = True)
      h3 = TextWrapTag("h3", p.title2, inline = True)
      divtag2.attachChild(h3, "child")
      divtag.attachChild(divtag2, "child")
      fixtag.attachChild(divtag, "tag-title2")      
      
      tag.attachChild(fixtag, "tag-title")
    else:
      divtag = SimpleTag ("div",[("class","field-container")])
      divtag2 = SimpleTag ("div",[("class","title")], inline = True)
      h3 = TextWrapTag("h3", p.title, inline = True)
      divtag2.attachChild(h3, "child")
      divtag.attachChild(divtag2, "child")
      tag.attachChild(divtag, "tag-title")
  
def _wrap_fieldcontainer(node, parent):
  nulltag = NullTag()
  #_add_headers(node, nulltag) #headers are ignored in ATTRACT web GUI
  _add_title(node, nulltag)  
  divtag = SimpleTag("div", [("class","field-container")], lines_before = 1 )
  _add_head(node, divtag)
  nulltag.attachChild(divtag, "child")
  parent.detachAndReplaceChild(node.name, nulltag)  
  divtag.attachChild(node, "child")
  return nulltag
  
def _container_tag(node, tag):
  assert node.space.htmlnode is None, (node.space.formpath, node.space.htmlnode)
  assert tag is not None, node.space.formpath
  Handler.passdown("html-make-concrete", node)
  node.space.htmlnode = tag  
  for childname in node.getChildren():
    child = node[childname]
    if child.space.complexity is not None:
      tip = child.space.htmltip
      if tip is None: tip = child.space.htmlnode      
      assert tip is not None, (child.space.formpath, child.parent.space.formpath, child.name, child.parent.name )
      #_add_headers(child, tip, begin=True) #headers are ignored in ATTRACT web GUI
      continue
    _wrap_fieldcontainer(child, node)
  
  node.morph("html-concrete-node")
  return Handler.CLAIM

def morph_formwrapper(node):
  assert node.space.cgi is not None
  attributes = [
    ("encType", "multipart/form-data"),
    ("action", node.space.cgi),
    ("method", "post"),
  ]
  if node.space.newtab:
    attributes.append(("target", "_blank"))
  comment = "<!-- END: form -->"
  formtag = SimpleTag("form", attributes, comment = comment)
  if node.space.hidden is not None: 
    for k,v in node.space.hidden.items():
      attributes = (
        ("type", "hidden"),
        ("name", k),
        ("value", str(v)),
      )
      hiddentag = SimpleTag("input", attributes, lines_before = 1)
      formtag.attachChild(hiddentag, "hiddenchild_" + k)
  node.space.htmlnode = None
  _container_tag(node, formtag)
  
def morph_level3(node):
  has_switch = False
  if node.space.group and node.space.group.elegroup:
    if hasattr(node.space.group.elegroup, "has_switch"):
      has_switch = node.space.group.elegroup.has_switch
  assert not has_switch, (node.name, node.parent.space.formpath)
  
  nulltag = NullTag()
  n = NullTag()
  _add_title(node, nulltag)
  nulltag.attachChild(n, "child")
  nulltag.space.htmltip = n
  return _container_tag(node, nulltag)

  
def _add_controls(node):
  
  controlchildren = []
  if not node.space.controltitle: return NullTag()
  t = ComplexTag (
    "div",
    [("class", "title")],
    [ TextWrapTag("h4", node.space.controltitle, inline = True) ],
    inline = True  
  )
  controlchildren.append(t)
  t = ComplexTag (
    "div",
    [("class", "button reload-icon reload-form-block")],
    [ComplexTag(
     "div",
     [("class","controls-tooltip")], 
     [TextWrapTag("p", "Restore default values",inline=True)], 
     inline=True
    )],    
    inline = True  
  )  
  controlchildren.append(t)
  if node.space.clone is not None and node.space.clone.button is not None:
    t = ComplexTag(
     "div", 
     [
      ("id", "remove-block-%s-0" % node.space.blockname),
      ("class", "button close-icon"),
     ],
     [ComplexTag(
      "div",
      [("class","controls-tooltip")], 
      [TextWrapTag("p", "Remove partner",inline=True)], 
      inline=True
     )],    
     inline = True,   
    )
    controlchildren.append(t)  
  return ComplexTag (
    "div",
    [("class", "controls")],
    controlchildren,
    lines_before = 1
  )
  
def morph_categorywrapper(node):
  space = node.space
  divtag = SimpleTag(
   "div", 
   [
     ("id", "page%d" % space.page),
     ("class", "category-container"),
   ], 
   lines_before = 1,
  )
  
  if space.html_description:
    desctexttag = StringTag(space.description)
  else:
    desctexttag = TextWrapTag("p", space.description, inline = True)
  texttag = ComplexTag (
    "div",
    [("class", "description")],    
    [
    TextWrapTag("h2", space.title, inline = True),
    desctexttag,
    ],     
  )
  desctag = ComplexTag (
   "div",
   [("class", "form-category-description")],
   [
     SimpleTag (
      "div",
      [("class", "large-icon %s" % space.icon)],
      inline = True
     ),
     texttag,
   ],     
   lines_after = 1,
  )
  if space.clones is not None and len(space.clones):
    children = []
    for clone in space.clones:
      if clone.button is not None:
	onClick = "cloneBlock('block-%s-',%d);" % (clone.name, clone.length)
	child = ComplexTag (
	  "div",
	  [("class", "field-container")],
	  [ComplexTag (
	    "div",
	    [("class", "field-item button align-right")],
	    [VoidTag(
	      "input",
	      [
		("type","button"),
		("value",clone.button),
		("id","block-%s-add" % clone.name ),
		("onClick", onClick),
	      ],
	      inline = True,
	    )],
	  )],
	)
	children.append(child)    
    controltag = ComplexTag (
      "div",
      [("class", "category-controls")],
      children,
      lines_before = 1,
      lines_after = 1,
    )
    desctag.attachChild(controltag, "control")
  divtag.attachChild(desctag, "description")
  
  commenttag = StringTag(
   "<!-- UI CONTAINER FOR FORM ELEMENTS OF A CATEGORY SELECTED IN CATEGORY NAVIGATION MENU -->",
  ) 
  divtag.attachChild(commenttag, "comment")
  
  contenttag = ComplexTag(
   "div",
   [("class","form-category")],
   comment = "<!-- END: UI container -->"
  )
  divtag.attachChild(contenttag, "content")
  divtag.space.htmltip = contenttag
  
  
  return _container_tag(node, divtag)
  

def morph_markup_fieldcontainer(node):      
  inside = False
  insides = []
  childnames0 = list(node.getChildren())
  for nr, childname in enumerate(node.getChildren()):
    space = node[childname].space
    fields = ("span",)
    p = _extract_fields(space, space.form, fields)
    if inside:
      if p.span: inside = None
      else: inside = False
    else:
      if not p.span and nr < len(node.getChildren()) - 1:
        inside = True
      if p.span: inside = None  
    insides.append(inside)
  insides2 = []
  for n in range(len(node.getChildren())):
    v = insides[n]
    if n < len(node.getChildren()) - 1 and v == True and insides[n+1] == None: v = False
    insides2.append(v)
  insides = insides2

  switch = node.space.switchtag
  if switch:
    switch.sendMessage("make-markup")
  tag = node.space.htmlnode
  node.space.htmlnode = None
  _container_tag(node, tag)
  node.morph("html-concrete-node")
  node.sendMessage("make-markup")
    
  childnames = list(node.getChildren())
  assert childnames == childnames0 #make-markup may not mess with children names
  state = False
  for inside, childname in zip(insides, childnames):
    if inside == True and state == False:
      divtag = SimpleTag("div", [("class","group-inline clearfix")] , 
       comment = " <!--END: group-inline -->", lines_inside_after = 1)
      child = node.detachAndReplaceChild(childname, divtag)
      divtag.attachChild(child, "left")
      state = True
    elif inside == False and state == True:  
      child = node.detachChild(childname)
      divtag.attachChild(child, "right")
      state = False
  
  return Handler.CLAIM  
  
def morph_level2(node):    
  Handler.passdown("html-make-concrete", node)
  divtag = SimpleTag("div", [("class", "level2")] , comment = " <!-- END: level2 -->", lines_before = 1)
  div2tag = SimpleTag("div", [("class", "group-container")], 
   comment = " <!-- END: group-container -->", lines_before = 1, lines_inside_before = 1)
  has_switch = False
  if node.space.group and node.space.group.elegroup:
    if hasattr(node.space.group.elegroup, "has_switch"):
      has_switch = node.space.group.elegroup.has_switch
  if has_switch:
    assert node.getChildren(), node.space.formpath
    firstchildname = node.getChildren()[0]    
    assert node[firstchildname].space.formtype == "switch", (node.space.formpath, firstchildname)
    firstchild = node.detachChild(firstchildname)
    divtag.attachChild(firstchild, "switch")    
    divtag.root = node.root
    firstchild.sendMessage("html-make-concrete")
    firstchild = _wrap_fieldcontainer(firstchild, divtag)
    node.space.switchtag = firstchild
  else:
    _add_title(node, div2tag)  
    nulltag = NullTag()
    div2tag.attachChild(nulltag, "tag-tip")
    div2tag.space.htmltip = nulltag
  divtag.attachChild(div2tag, "container")  
  divtag.space.htmltip = div2tag
  node.space.htmlnode = divtag
  node.morph("html-concrete-fieldcontainer")

def morph_cloneblockwrapper(node):  
  return morph_level1_blockwrapper(node)

def morph_clonewrapper(node):
  return _container_tag(node, NullTag())
  
def morph_level1_blockwrapper(node):  
  space = node.space

  controltag = _add_controls(node)
  containertag = SimpleTag (
   "div",
   [("class","level1-container")],
   lines_before = 1,
   comment = "<!-- END: level1-container -->"
  )
  
  level1tag = ComplexTag (
   "div",
   [
    ("id","block-%s-0" % node.space.blockname),
    ("class", "level1 active")
   ],
   [ controltag, containertag ],
   lines_before = 1,
   lines_after = 1,
   comment = "<!-- END: level1 -->"
  )
  level1tag.space.htmltip = containertag
  Handler.passdown("html-make-concrete", node)
  return _container_tag(node, level1tag)
    
