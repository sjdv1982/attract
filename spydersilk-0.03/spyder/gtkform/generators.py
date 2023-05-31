from functools import partial 

def generator_textarea(name, id, objs, forms, wrappings, defaultvalue):
  form = forms[0]
  expand = True
  if hasattr(form, "expand"): expand = form.expand
  if hasattr(form, "name"): name = form.name
  
  ret = []
  ret.append('<object class="GtkLabel" id="_head-%s">' % id)
  ret.append('  <property name="visible">True</property>')
  ret.append('  <property name="hexpand">%s</property>' % expand)
  ret.append('  <property name="vexpand">%s</property>' % expand)
  ret.append('  <property name="label" translatable="yes">%s</property>' % name)  
  ret.append('  <property name="can_focus">False</property>')
  ret.append('</object>')
  ret.append('<packing>')
  ret.append('  <property name="expand">%s</property>' % expand)
  ret.append('  <property name="fill">%s</property>' % expand)        
  ret.append('  <property name="left_attach">0</property>')
  ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
  ret.append('  <property name="width">2</property>')
  ret.append('  <property name="height">1</property>')            
  ret.append('</packing>')
  ret.append('</child><child>')
  ret.append('<object class="GtkTextView" id="_widget-%s">' % id)
  ret.append('  <property name="visible">True</property>')
  ret.append('  <property name="hexpand">%s</property>' % expand)
  ret.append('  <property name="vexpand">%s</property>' % expand)
  ret.append('  <property name="can_focus">True</property>')
  ret.append('  <property name="buffer">_textbuffer-%s</property>' % id)  
  ret.append('</object>')
  ret.append('<packing>')
  ret.append('  <property name="expand">%s</property>' % expand)
  ret.append('  <property name="fill">%s</property>' % expand)        
  ret.append('  <property name="left_attach">0</property>')
  ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
  ret.append('  <property name="width">2</property>')
  ret.append('  <property name="height">1</property>')            
  ret.append('</packing>')

  ret2 = []
  if defaultvalue is None: defaultvalue = ""
  ret2.append('<object class="GtkTextBuffer" id="_textbuffer-%s">' % id)
  ret2.append('  <property name="text">%s</property>' % defaultvalue)  
  ret2.append('</object>')  
  return ret, ret2

def generator_expander(name, id, objs, forms, wrappings, mode):
  if mode == "start": 
    ret = []
    ret.append('<object class="GtkExpander" id="_exp-%s">' % id)
    ret.append('  <property name="visible">True</property>')
    ret.append('  <property name="can_focus">True</property>')
    return ret, None
  else: 
    ret = []
    ret.append("</object> <!-- GtkExpander -->")
    ret.append("<packing>")
    form = forms[0]
    expand = True
    if hasattr(form, "expand"): expand = form.expand
    ret.append('  <property name="expand">%s</property>' % expand)
    ret.append('  <property name="fill">%s</property>' % expand)        
    ret.append('  <property name="left_attach">0</property>')
    ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
    ret.append('  <property name="width">2</property>')
    ret.append('  <property name="height">1</property>')            
    ret.append("</packing>")
    return ret, None

header_count = {}

def generator_header(name, id, objs, forms, wrappings, mode, args):
  if mode != "start": return None, None
  if id not in header_count: header_count[id] = 0
  header_count[id] += 1
  idstr = str(id) + "-" + str(header_count[id])
  
  form = forms[0]
  expand = True
  if hasattr(form, "expand"): expand = form.expand
  
  ret = []
  ret.append('<child> <!-- generate_header GtkLabel -->')
  ret.append('  <object class="GtkLabel" id="_header-%s">' % idstr)
  ret.append('    <property name="visible">True</property>')
  ret.append('    <property name="hexpand">%s</property>' % expand)
  ret.append('    <property name="vexpand">%s</property>' % expand)
  ret.append('    <property name="can_focus">False</property>')
  ret.append('    <property name="label" translatable="yes">%s</property>' % args[0]) 
  #ret.append('    <property name="use_markup">True</property>') #TODO  
  ret.append("  </object> <!-- generate_header GtkLabel -->")
  ret.append("  <packing>")
  ret.append('    <property name="expand">%s</property>' % expand)
  ret.append('    <property name="fill">%s</property>' % expand)        
  ret.append('    <property name="left_attach">0</property>')
  ret.append('    <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
  ret.append('    <property name="width">2</property>')
  ret.append('    <property name="height">1</property>')            
  ret.append("  </packing>")
  ret.append('</child> <!-- generate_header GtkLabel -->')

  return ret, None

button_count = {}

def generator_button(name, id, objs, forms, wrappings, mode, args):
  assert len(args) == 2, len(args)
  assert args[-1] in ("after", "before")
  if args[-1] == "before":
    return generator_pre_button(name, id, objs, forms, wrappings, mode, args[0])
  elif args[-1] == "after":
    return generator_post_button(name, id, objs, forms, wrappings, mode, args[0])

def generator_pre_button(name, id, objs, forms, wrappings, mode, text):
  if mode != "start": return None, None
  if id not in button_count: button_count[id] = 0
  button_count[id] += 1
  idstr = str(id) + "-" + str(button_count[id])
  
  form = forms[0]
  expand = True
  if hasattr(form, "expand"): expand = form.expand
  
  ret = []
  ret.append('<child>')
  ret.append('  <object class="GtkButton" id="_button-%s">' % idstr)
  ret.append('    <property name="visible">True</property>')
  ret.append('    <property name="hexpand">%s</property>' % expand)
  ret.append('    <property name="vexpand">%s</property>' % expand)
  ret.append('    <property name="can_focus">False</property>')
  ret.append('    <property name="label" translatable="yes">%s</property>' % text) 
  #ret.append('    <property name="use_markup">True</property>') #TODO  
  ret.append("  </object>")
  ret.append("  <packing>")
  ret.append('    <property name="expand">%s</property>' % expand)
  ret.append('    <property name="fill">%s</property>' % expand)        
  ret.append('    <property name="left_attach">0</property>')
  ret.append('    <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
  ret.append('    <property name="width">2</property>')
  ret.append('    <property name="height">1</property>')            
  ret.append('  </packing>')
  ret.append('</child>')

  return ret, None

def generator_post_button(name, id, objs, forms, wrappings, mode,text):
  if mode == "start": return None, None
  if id not in button_count: button_count[id] = 0
  button_count[id] += 1
  idstr = str(id) + "-" + str(button_count[id])
  
  form = forms[0]
  expand = True
  if hasattr(form, "expand"): expand = form.expand
  
  ret = []
  ret.append('<child>')
  ret.append('  <object class="GtkButton" id="_button-%s">' % idstr)
  ret.append('    <property name="visible">True</property>')
  ret.append('    <property name="hexpand">%s</property>' % expand)
  ret.append('    <property name="vexpand">%s</property>' % expand)
  ret.append('    <property name="can_focus">False</property>')
  ret.append('    <property name="label" translatable="yes">%s</property>' % text) 
  #ret.append('    <property name="use_markup">True</property>') #TODO  
  ret.append("  </object>")
  ret.append("  <packing>")
  ret.append('    <property name="expand">%s</property>' % expand)
  ret.append('    <property name="fill">%s</property>' % expand)        
  ret.append('    <property name="left_attach">0</property>')
  ret.append('    <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
  ret.append('    <property name="width">2</property>')
  ret.append('    <property name="height">1</property>')            
  ret.append('  </packing>')
  ret.append('</child>')

  return ret, None
  
def generator_box(ori, name, id, objs, forms, wrappings, mode):
  if mode == "start": 
    ret = []
    ret.append('<object class="GtkBox" id="_box-%s">' % id)
    ret.append('  <property name="visible">True</property>')
    ret.append('  <property name="orientation">%s</property>' % ori)
    return ret, None
  else: 
    ret = []
    ret.append("</object> <!-- GtkBox %s-->" % id)
    ret.append("<packing>")
    form = forms[0]
    expand = True
    if hasattr(form, "expand"): expand = form.expand
    ret.append('  <property name="expand">%s</property>' % expand)
    ret.append('  <property name="fill">%s</property>' % expand)        
    ret.append('  <property name="left_attach">0</property>')
    ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
    ret.append('  <property name="width">2</property>')
    ret.append('  <property name="height">1</property>')            
    ret.append("</packing>")
    return ret, None

def generator_grid_toplevel(name, id, objs, forms, wrappings, mode):
  if wrappings == 0: return generator_grid(name, id, objs, forms, wrappings, mode)
  else: return None, None

def generator_grid(name, id, objs, forms, wrappings, mode):
  if mode == "start": 
    ret = []
    form = forms[0]
    expand = True
    if hasattr(form, "expand"): expand = form.expand    
    ret.append('<object class="GtkGrid" id="_box-%s">' % id)
    ret.append('  <property name="visible">True</property>')
    ret.append('  <property name="hexpand">%s</property>' % expand)
    ret.append('  <property name="vexpand">%s</property>' % expand)
    ret.append('  <property name="columns">2</property>')
    ret.append('  <property name="n_columns">2</property>')
    return ret, None
  else: 
    ret = []
    ret.append("</object> <!-- GtkGrid %s -->" % id)
    ret.append("<packing>")
    form = forms[0]
    expand = True
    if hasattr(form, "expand"): expand = form.expand
    ret.append('  <property name="expand">%s</property>' % expand)
    ret.append('  <property name="fill">%s</property>' % expand)        
    ret.append('  <property name="left_attach">0</property>')
    ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
    ret.append('  <property name="width">2</property>')
    ret.append('  <property name="height">1</property>')            
    ret.append("</packing>")
    return ret, None

def generator_grid2(name, id, objs, forms, wrappings, mode):
  if mode == "start": 
    ret = []
    ret.append('<object class="GtkGrid" id="_box-%s">' % id)
    ret.append('  <property name="visible">True</property>')
    ret.append('  <property name="columns">2</property>')
    ret.append('  <property name="n_columns">2</property>')
    return ret, None
  else: 
    ret = []
    ret.append("</object> <!-- GtkGrid %s -->" % id)
    ret.append("<packing>")
    form = forms[0]
    expand = True
    if hasattr(form, "expand"): expand = form.expand
    ret.append('  <property name="expand">%s</property>' % expand)
    ret.append('  <property name="fill">%s</property>' % expand)        
    ret.append('  <property name="left_attach">0</property>')
    ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
    ret.append('  <property name="width">2</property>')
    ret.append('  <property name="height">1</property>')            
    ret.append("</packing>")
    return ret, None

def generator_child(name, id, objs, forms, wrappings, mode):
  if mode == "start": 
    return "<child>", None
  else: 
    return "</child> <!-- %s -->" % id, None

def generator_wrap_child(gen, *args):
  mode = args[-1]
  id = args[1]
  wrappings = args[-2]
  forms = args[-3]
  assert mode in ("start", "end")
  if wrappings == 0: #toplevel form, no outerwrap or head elements
    return gen(*args)
  ret = []
  if mode == "start" :ret.append("<child>")
  r = gen(*args)[0]
  for rr in r:
    ret.append("  " + rr)
  if mode == "end" :ret.append("</child> <!-- %s -->" % id)
  return ret, None

def wrap_child(gen):
  return partial(generator_wrap_child, gen)

def generator_head(childtype, name, id, objs, forms, wrappings, name2):
  ret = []
  if childtype:
    ret.append('<child type="label">')
  else:
    ret.append('<child>')
  ret.append('  <object class="GtkLabel" id="_head-%s">' % id)
  ret.append('    <property name="visible">True</property>')
  ret.append('    <property name="label" translatable="yes">%s</property>' % name2)
  ret.append('  </object>')
  ret.append('  <packing>')
  form = forms[0]
  expand = True
  if hasattr(form, "expand"): expand = form.expand
  ret.append('    <property name="expand">%s</property>' % expand)
  ret.append('  <property name="fill">%s</property>' % expand)        
  if not childtype:
    ret.append('  <property name="left_attach">0</property>')
    ret.append('  <property name="top_attach">!!NEW-MEMBERPOS!!</property>')
    ret.append('  <property name="width">2</property>')
    ret.append('  <property name="height">1</property>')            
  ret.append('  </packing>')
  ret.append('</child>')
  return ret, None

def generator_column_or_grid(name, id, objs, forms, wrappings, mode):
  form = forms[0]
  if hasattr(form, "subtype") and form.subtype == "column":
    return generator_box("horizontal", name, id, objs, forms, wrappings, mode)
  else:
    return generator_grid2(name, id, objs, forms, wrappings, mode)      

def generator_group_expander(name, id, objs, forms, wrappings, group, groupindex, mode):
  if mode in ("pre", "post"): return None
  if id is None: groupid = "group" + str(groupindex)
  else: groupid = id + "-group" + str(groupindex)
  if mode == "start":
    ret = ["<child>"]
    ret1 = generator_expander(groupid, groupid, objs, forms, wrappings, mode)
    ret2 = generator_head(True, groupid, groupid, objs, forms, wrappings, group)
    ret3 = generator_grid(groupid, groupid, objs, forms, wrappings, mode)
    return (ret + ret1[0] + ret2[0] + ret + ret3[0]),None
  elif mode == "end":
    ret = ["</child>"]
    ret1 = generator_expander(groupid, groupid, objs, forms, wrappings, mode)
    ret3 = generator_grid(groupid, groupid, objs, forms, wrappings, mode)
    return (ret3[0] + ret + ret1[0] + ret), None
    

def init():
  from .xml import set_generator, set_element_generator, set_group_generator, set_phase_generator

  set_element_generator("textarea", generator_textarea)
  
  set_generator("outerwrap", 0, None)
  set_phase_generator("header", 0, "pre", generator_header)    
  set_phase_generator("button", 0, "pre", generator_button) 
  set_generator("group", 0, None) #never used
  set_generator("midwrap", 0, generator_child)
  set_generator("head", 0, None) #head is already part of generated element
  set_generator("innerwrap", 0, None)

  set_generator("outerwrap", 1, None)
  set_phase_generator("header", 1, "pre", generator_header)    
  set_phase_generator("button", 1, "inner", generator_button) 
  set_generator("group", 1, None) #never used
  set_generator("midwrap", 1, wrap_child(generator_grid))  
  set_generator("head", 1, partial(generator_head, False))
  set_generator("innerwrap", 1, generator_grid_toplevel)
  
  set_generator("outerwrap", 2, None)
  set_phase_generator("header", 2, "pre", generator_header) 
  set_phase_generator("button", 2, "inner", generator_button)      
  set_generator("group", 2, generator_group_expander)
  set_generator("midwrap", 2, wrap_child(generator_expander))
  set_generator("head", 2, partial(generator_head, True))
  set_generator("innerwrap", 2, wrap_child(generator_column_or_grid))
  
  for n in range(3, 20):
    set_generator("outerwrap", n, None)
    set_phase_generator("header", n, "pre", generator_header)        
    set_phase_generator("button", n, "inner", generator_button)      
    set_generator("group", n, generator_group_expander)
    set_generator("midwrap", n, wrap_child(generator_expander))
    set_generator("head", n, partial(generator_head, True))
    set_generator("innerwrap", n, wrap_child(generator_grid))
  
"""
  set_generator("button", 0, generator_button_0)  
  
  set_generator("button", 1, generator_button_1)  

  set_generator("button", 2, generator_button_2)  

  for n in range(3, 20):
    set_generator("button", n, generator_button_2)  

def generator_button_0(*args):
  phase = args[-2]
  if phase == "pre": return generator_pre_button(*args)
  elif phase == "outer": return generator_post_button(*args)
  else: raise Exception(phase)

def generator_button_1(*args):
  phase = args[-2]
  if phase in ("inner", "pre"): return generator_pre_button(*args)
  elif phase == "post": return generator_post_button(*args)
  else: raise Exception(phase)

def generator_button_2(*args):
  place = args[-1][1]
  phase = args[-2]
  if place == "before": return generator_pre_button(*args)
  elif place == "after": return generator_post_button(*args)
  else: raise Exception(place)

"""  
