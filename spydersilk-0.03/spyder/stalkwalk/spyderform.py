from .Node import Root
from . import Handler
from .Space import Space
from .Handler import Constructor, CMethod, XMethod, Method, UniMethod, PassDown, NullHandler
import spyder, Spyder
from Spyder import Object, File, String, Integer, Bool, Float

def init_spyderform(node):
  params = node.parent.space.params
  space = node.space
  
  form = params.spyderform
  obj = params.obj
  space.obj = obj
  
  if hasattr(form, "value"):
    defaultvalue = form.value
  elif obj is not None:
    defaultvalue = obj
  elif hasattr(form, "default"):
    defaultvalue = form.default
  else:
    defaultvalue = None
  space.defaultvalue = defaultvalue
  
  curr_membernames = None
  curr_members = None
  curr_objs = None
  if form.arraycount > 0:
    is_array = True
    is_none = False   
    if hasattr(form, "type"):
      if form.type not in ("none",None): is_array = False
      if form.type == "none": is_none = True
    if is_array:
      if is_none is False:
        n0, typnam = node.name, None
        if hasattr(form, "_typetree"):
          typnam = form._typetree.typename + "Array" * form.arraycount
        if not hasattr(form, "length"):
          raise TypeError("Cannot generate for \"%s\"(%s): member is an Array, but no length has been specified" % (n0, typnam))      
        if not isinstance(form.length, int) or form.length <= 0:
          raise TypeError("Cannot generate for \"%s\"(%s): member is an Array, but its length is not a positive integer" % (n0, typnam))                  
        curr_membernames = []
        curr_members = []
        curr_objs = []
        for nr in range(form.length):
          curr_membernames.append(nr)
          curr_members.append(form[nr])
          mobj = None
          if defaultvalue is not None and len(defaultvalue) > nr:
            mobj = defaultvalue[nr]
          curr_objs.append(mobj)  
  elif form.get_membernames():
    curr_membernames = form.get_membernames()
    curr_members = []
    curr_objs = []
    for membername in curr_membernames:
      curr_members.append(getattr(form, membername, None))
      mobj = None
      if defaultvalue is not None: 
        mobj = getattr(defaultvalue, membername, None)
      curr_objs.append(mobj)

  space.curr_membernames = curr_membernames
  space.curr_members = curr_members
  space.curr_objs = curr_objs
  attrs = (
   "memberorder", "group", "options", "optiontitles"
  )
  for a in attrs:
    if hasattr(form, a):
      setattr(space,a, getattr(form,a)) 
  if hasattr(form, "group"):
    space.group = form.group     
  space.elegroups = form._groups  
  space.form = Space(form._props.copy())
  if "_othertokens" in form._props:
    space.othertokens = list(form._props["_othertokens"])
  if hasattr(form, "type"):
    if form.type is None:
      space.formtype = "none"
    else:
      space.formtype = form.type 
  if hasattr(form, "subtype"):
    space.formsubtype = form.subtype 
  if hasattr(form, "typename"):
    space.formtypename = form.typename 
  space.formname = node.name
  if hasattr(form, "name"):
    space.formname = form.name
    
  if curr_membernames is not None:    
    if hasattr(form, "memberorder"):
      mems = []
      for mnr, m in enumerate(form.memberorder):
        assert m not in form.memberorder[mnr+1:]
        if m not in curr_membernames: continue
        mems.append(curr_membernames.index(m))
      for mnr, membername in enumerate(curr_membernames):        
        if membername not in form.memberorder:
          mems.append(mnr)
      curr_membernames = [curr_membernames[i] for i in mems]      
      curr_members = [curr_members[i] for i in mems]
      curr_objs = [curr_objs[i] for i in mems]

  if curr_membernames is not None:          
    for membername, mform, mobj in zip(curr_membernames, curr_members, curr_objs):
      node.space.params.spyderform = mform
      node.space.params.obj = mobj
      node.addChild("spyderform", membername)
      node.space.params.clear()
  
def determine_leaves(node):
  space = node.space
  if node.getChildren(): 
    if space.formtype is None:
      return Handler.NO_MATCH
  
  typ = Spyder.String
  if space.formtypename is not None:
    try:
      typ = getattr(Spyder, space.formtypename)
    except AttributeError:
      pass
  
  formtype = space.formtype
  if formtype is None:
    formtype = "text"
    if issubclass(typ, File):
      formtype = "file"
    elif issubclass(typ, Integer) or issubclass(typ, Float):
      formtype = "spin"
      
    elif issubclass(typ, Bool):
      formtype = "checkbox"
      
    if space.options is not None:
      assert formtype != "file"
      formtype = "option"
      options = space.options      
      optiontitles = options
      if space.optiontitles is not None:
        optiontitles = space.optiontitles
        assert len(options) == len(optiontitles), (node.name, len(options), len(optiontitles))
      space.optiontitles = optiontitles
      
    space.formtype = formtype
  
  formsubtype = space.formsubtype
  if formsubtype is None:    
    if issubclass(typ, Integer):
      formsubtype = "int"
    elif issubclass(typ, Float):
      formsubtype = "float"
    space.formsubtype = formsubtype
    
  for childname in list(node.getChildren()):
    node.removeChild(childname)
  
  #print("LEAF!", node.name, space.formtype, space.formsubtype)
  node.morph("spyderwidget")
  return Handler.CLAIM
  
class SpyderformRoot(Root): pass 

def determine_wrapping(node):
  space = node.space
  if node.parent.space.wrapping is None:
    space.wrapping = 0
  else:
    space.wrapping = node.parent.space.wrapping + 1
  for child in node.getChildren():
    node[child].sendMessage("determine-wrapping")
  return Handler.CLAIM
  
def wrapping_bounce(node):
  space = node.space
  space.wrapping = node.parent.space.wrapping + 1
  node.parent.sendMessage("update-complexity", args=(1,))
  return Handler.CLAIM

def update_complexity(node, complexity):
  space = node.space
  if space.complexity is None: space.complexity = 0
  if complexity > space.complexity:
    space.complexity = complexity
    node.parent.sendMessage("update-complexity", args=(complexity+1,))
  return Handler.CLAIM

def pass_formpath(node, currpath):
  if node.name is not None:
    currpath = currpath + (node.name,)
  node.space.formpath = currpath
  for child in node.getChildren():
    node[child].sendMessage("set-formpath", args = (currpath,) )
  return Handler.CLAIM

def set_formpath(node, currpath):
  node.space.formpath = currpath + (node.name,)
  return Handler.CLAIM
  
def prune_tree(node):
  for child in list(node.getChildren()):
    node[child].sendMessage("prune-tree")  
  if not node.getChildren():
    node.parent.removeChild(node.name)
SpyderformRoot.addHandler( NullHandler("root", "update-complexity" ) )

SpyderformRoot.addHandler( Constructor("spyderform", init_spyderform) )
SpyderformRoot.addHandler( PassDown("spyderform", "determine-leaves") )
SpyderformRoot.addHandler( CMethod("spyderform", "determine-leaves", determine_leaves ) )
from .formgroup import remove_member, determine_groups, prune_groups, make_groups, init_group, complete_group
SpyderformRoot.addHandler( XMethod("spyderform", "grouping-remove-member", remove_member ) )
SpyderformRoot.addHandler( XMethod("spyderform", "determine-groups", determine_groups ) )
SpyderformRoot.addHandler( XMethod("spyderform", "prune-groups", prune_groups ) )
SpyderformRoot.addHandler( Method("spyderform", "determine-wrapping", determine_wrapping ) )
SpyderformRoot.addHandler( Method("spyderform", "update-complexity", update_complexity ) )
import functools
grouphandler = functools.partial(make_groups, "spydergroup")
SpyderformRoot.addHandler( XMethod("spyderform", "make-groups", grouphandler ) )
SpyderformRoot.addHandler( PassDown("spyderform", "complete-groups") )
SpyderformRoot.addHandler( Method("spyderform", "prune-tree", prune_tree ) )
SpyderformRoot.addHandler( Method("spyderform", "set-formpath", pass_formpath ) )

SpyderformRoot.addHandler( Constructor("spydergroup", init_group) )
SpyderformRoot.addHandler( NullHandler("spydergroup", "make-groups" ) )
SpyderformRoot.addHandler( XMethod("spydergroup", "complete-groups", complete_group ) )
SpyderformRoot.addHandler( Method("spydergroup", "prune-tree", prune_tree ) )
SpyderformRoot.addHandler( Method("spydergroup", "determine-wrapping", determine_wrapping ) )
SpyderformRoot.addHandler( Method("spydergroup", "update-complexity", update_complexity ) )

SpyderformRoot.addHandler( NullHandler("spyderwidget", "determine-groups" ) )
SpyderformRoot.addHandler( NullHandler("spyderwidget", "prune-groups" ) )
SpyderformRoot.addHandler( NullHandler("spyderwidget", "make-groups" ) )
SpyderformRoot.addHandler( NullHandler("spyderwidget", "complete-groups" ) )
SpyderformRoot.addHandler( NullHandler("spyderwidget", "prune-tree" ) )
SpyderformRoot.addHandler( Method("spyderwidget", "determine-wrapping", wrapping_bounce ) )
SpyderformRoot.addHandler( Method("spyderwidget", "set-formpath", set_formpath ) )

def build(root, spyderform, obj = None):
  root.space.params.spyderform = spyderform
  root.space.params.obj = obj  
  root.addChild("spyderform", "top")
  root.space.params.clear()
  root["top"].sendMessage("determine-leaves")
  root["top"].sendMessage("set-formpath", args = ((),))
  root["top"].sendMessage("determine-groups")
  root["top"].sendMessage("prune-groups")
  root["top"].sendMessage("make-groups") 
  root["top"].sendMessage("complete-groups") 
  root["top"].sendMessage("prune-tree") 
  root["top"].sendMessage("determine-wrapping")
  