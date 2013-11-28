from .concrete import _extract_fields
from spyder.stalkwalk import Handler, Space
from spyder.stalkwalk.hamlet.tag import SimpleTag, VoidTag, ComplexTag, ComplexVoidTag, TextWrapTag, StringTag, NullTag

  
def _leaf_tag(node, tag):
  node.space.htmlnode = tag
  Handler.passdown("html-make-concrete", node)  
  node.morph("html-concrete-node")
  return Handler.CLAIM

def _div_leaf_tag(node,tag, divclass):
  divtag = SimpleTag("div", [("class", divclass)])
  divtag.attachChild(tag, "child")
  node.space.htmlnode = divtag
  Handler.passdown("html-make-concrete", node)  
  node.morph("html-concrete-node")
  return Handler.CLAIM
  
def morph_inputtag_switch(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "switch": return Handler.NO_MATCH
  fields = ("defaultvalue",)
  p = _extract_fields(space, space.form, fields)
  checked = None
  attributes = [
    ("type", "checkbox"), 
    ("name", fp),
    ("id", fp),
    ("value", ""),
  ]
  if p.defaultvalue: attributes.append(("checked","checked"))
  inputtag = VoidTag("input", attributes)
  divtag = ComplexTag (
    "div",
    [ ("class", "field-item switch align-left"), ],
    ( 
     inputtag,
     TextWrapTag("label", "Off", [ ("for", fp) ], inline = True ),
     TextWrapTag("label", "On",  [ ("for", fp) ], inline = True ),
     SimpleTag("span", [ ("class", "toggle") ], inline = True)
    ),  
  )
  return _leaf_tag(node, divtag)
  
def morph_inputtag_text(node):
  space = node.space
  if space.formtype not in ("text", "password"): return Handler.NO_MATCH
  return _morph_inputtag_text(node, space.formtype)

def _morph_inputtag_text(node, typ):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  fields = ("defaultvalue", "length", "placeholder")
  p = _extract_fields(space, space.form, fields)
  if p.defaultvalue is None: p.defaultvalue = ""
  p.defaultvalue = str(p.defaultvalue) 
  if p.length is not None: p.length = str(p.length)
  attributes = [
    ("type", typ), 
    ("name", fp),
    ("id", fp),
    ("value", p.defaultvalue),
    ("maxlength", p.length),
  ]  
  attributes.append(("placeholder", p.placeholder))
  inputtag = VoidTag("input", attributes)
  return _div_leaf_tag(node, inputtag, "field-item %s align-right" % typ)
  
def morph_textarea(node):
  space = node.space
  if space.formtype != "textarea": return Handler.NO_MATCH
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  fields = ("defaultvalue", "rows")
  p = _extract_fields(space, space.form, fields)
  if p.defaultvalue is None: p.defaultvalue = ""
  p.defaultvalue = str(p.defaultvalue) 
  if p.rows is not None: p.rows = str(p.rows)
  attributes = [
    ("name", fp),
    ("rows", p.rows),
    ("id", fp),
  ]  
  textareatag = TextWrapTag("textarea", p.defaultvalue, attributes, inline=True)
  return _div_leaf_tag(node, textareatag, "field-item textarea align-right")
  
def morph_inputtag_number(node):
  space = node.space
  if space.formtype not in ("number", "spin"): return Handler.NO_MATCH
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  fields = ("defaultvalue", "min", "max", "step", "digits", "placeholder")
      
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("type", "number"), 
    ("name", fp),
    ("id", fp),    
  ]
  if p.defaultvalue is not None:
    attributes.append (
     ("value", str(p.defaultvalue), {"space_before_assign": 1, "space_after_assign": 1} )
    )  
  
  if p.min is not None:
    if p.digits is not None:
      formatstr = "%." + str(p.digits) + "f"
      minstr = formatstr % p.min
    elif isinstance(p.min, float):
      minstr = "%f" % p.min
    else:
      minstr = "%d" % p.min
    attributes.append(("min", minstr))  
  if p.max is not None:
    if p.digits is not None:
      formatstr = "%." + str(p.digits) + "f"
      maxstr = formatstr % p.max
    elif isinstance(p.max, float):
      maxstr = "%f" % p.max
    else:
      maxstr = "%d" % p.max
    attributes.append(("max", maxstr))
  if p.step is not None:
    if p.digits is not None:
      formatstr = "%." + str(p.digits) + "f"
      stepstr = formatstr % p.step
    elif isinstance(p.step, float):
      stepstr = "%f" % p.step
    else:
      stepstr = "%d" % p.step
    attributes.append(("step", stepstr))  
  
  attributes.append(("placeholder", p.placeholder))
  
  inputtag = VoidTag("input", attributes)
  return _div_leaf_tag(node, inputtag, "field-item number align-right")
  
def morph_selecttag(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "option": return Handler.NO_MATCH  
  fields = ("defaultvalue",)
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("name", fp),
    ("id", fp),
  ]  
  options, optiontitles = space.options, space.optiontitles
  assert len(options) and len(options) == len(optiontitles)
  children = []
  if p.defaultvalue is None:
    child = SimpleTag("option", (("value",""),"selected"), inline=True )
    children.append(child)
  for opt, optitle in zip(options, optiontitles):  
    attrs = (("value",str(opt)),)
    if opt == p.defaultvalue: attrs = (("value",str(opt)),"selected") 
    child = TextWrapTag("option", optitle, attrs, inline=True)
    children.append(child)
  inputtag = ComplexTag("select", attributes, children)
  
  return _div_leaf_tag(node, inputtag, "field-item select align-right")

def morph_inputtag_file(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "file": return Handler.NO_MATCH
  fields = ("placeholder","defaultvalue")
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("type", "file"), 
    ("name", fp),
    ("id", fp),
    ("placeholder", p.placeholder),
  ]  
  inputtag = VoidTag("input", attributes)
  if p.defaultvalue is not None:
    space.formname = "<em><mark>&lt;Embedded resource&gt;:</mark>&nbsp;&nbsp;&nbsp;</em>" + space.formname
  return _div_leaf_tag(node, inputtag, "field-item file align-right")

def morph_inputtag_checkbox(node):
  space = node.space
  assert space.formpath is not None
  fp = "-".join([str(s) for s in space.formpath[1:]])
  if space.formtype != "checkbox": return Handler.NO_MATCH
  fields = ("defaultvalue",)
  p = _extract_fields(space, space.form, fields)
  attributes = [
    ("type", "checkbox"), 
    ("name", fp),
    ("id", fp),
    ("value", ""),
  ]
  if p.defaultvalue: attributes.append(("checked","checked"))
  inputtag = VoidTag("input", attributes)
  return _div_leaf_tag(node, inputtag, "field-item checkbox align-right")
  
def morph_inputtag_unknown(node):
  raise ValueError("Unknown form type", (node.space.formtype, node.space.formpath))