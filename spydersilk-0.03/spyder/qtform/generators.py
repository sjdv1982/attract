from functools import partial 

def generator_textarea(name, id, objs, forms, wrappings, defaultvalue):
  ret = []
  form = forms[0]
  if hasattr(form, "name"): name = form.name

  ret.append('<item row="!!NEW-MEMBERPOS!!" column="0">')
  ret.append(' <widget class="QLabel" name="_head-%s">' % id)
  ret.append('  <property name="text">')
  ret.append('   <string>%s</string>' % name) 
  ret.append('  </property>')
  ret.append(' </widget>')
  ret.append('</item>')
  ret.append('<item row="!!MEMBERPOS!!" column="1">')
  ret.append(' <widget class="QPlainTextEdit" name="_widget-%s">' % id)
  ret.append('  <property name="sizePolicy">')  
  ret.append('   <sizepolicy hsizetype="Expanding" vsizetype="Expanding">')
  ret.append('    <horstretch>0</horstretch>')
  ret.append('    <verstretch>0</verstretch>')
  ret.append('   </sizepolicy>')   
  ret.append('  </property>')  
  ret.append('  <property name="minimumSize">')
  ret.append('   <size>')
  ret.append('    <width>600</width>')
  ret.append('    <height>0</height>')
  ret.append('   </size>')
  ret.append('  </property>')
  
  if defaultvalue is not None:
    ret.append('  <property name="plainText">')
    ret.append('   <string>%s</string>' % defaultvalue) 
    ret.append('  </property>')
  ret.append(' </widget>')
  ret.append('</item>')

  return ret, [], [], []

def generator_expander(name, id, objs, forms, wrappings, mode):
  if wrappings == 0: return None, None, None, None
  ret = []
  con = []
  if mode == "start": 
    ret.append('<layout class="QVBoxLayout" name="_exp-%s">' % id)
    #ret.append(' <item row="!!NEW-MEMBERPOS!!" column="0" colspan="2">')
    ret.append(' <item>')
    ret.append('  <layout class="QHBoxLayout" name="_exph-%s">' % id)
    ret.append('   <item>')
    ret.append('    <widget class="QPushButton" name="_exp_bplus-%s">' % id)
    ret.append('     <property name="text">')
    ret.append('      <string>+</string>')
    ret.append('     </property>')
    ret.append('    </widget>')
    ret.append('   </item>')
    ret.append('   <item>')
    ret.append('    <widget class="QPushButton" name="_exp_bminus-%s">' % id)
    ret.append('     <property name="text">')
    ret.append('      <string>-</string>')
    ret.append('     </property>')
    ret.append('     <property name="__hidden_at_startup__" stdset="0">')
    ret.append('      <bool>true</bool>')
    ret.append('     </property>')
    ret.append('    </widget>')
    ret.append('   </item>')
    ret.append('   <item>')
    ret.append('<!-- expansion -->')

    con.append(' <connection>')
    con.append('  <sender>_exp_bplus-%s</sender>' % id)
    con.append('  <signal>clicked()</signal>')
    con.append('  <receiver>_exp_bplus-%s</receiver>' % id)
    con.append('  <slot>hide()</slot>')    
    con.append(' </connection>')
    con.append(' <connection>')
    con.append('  <sender>_exp_bplus-%s</sender>' % id)
    con.append('  <signal>clicked()</signal>')
    con.append('  <receiver>_exp_bminus-%s</receiver>' % id)
    con.append('  <slot>show()</slot>')    
    con.append(' </connection>')
    con.append(' <connection>')
    con.append('  <sender>_exp_bplus-%s</sender>' % id)
    con.append('  <signal>clicked()</signal>')
    con.append('  <receiver>_expw-%s</receiver>' % id)
    con.append('  <slot>show()</slot>')    
    con.append(' </connection>')

    con.append(' <connection>')
    con.append('  <sender>_exp_bminus-%s</sender>' % id)
    con.append('  <signal>clicked()</signal>')
    con.append('  <receiver>_exp_bminus-%s</receiver>' % id)
    con.append('  <slot>hide()</slot>')    
    con.append(' </connection>')
    con.append(' <connection>')
    con.append('  <sender>_exp_bminus-%s</sender>' % id)
    con.append('  <signal>clicked()</signal>')
    con.append('  <receiver>_exp_bplus-%s</receiver>' % id)
    con.append('  <slot>show()</slot>')    
    con.append(' </connection>')
    con.append(' <connection>')
    con.append('  <sender>_exp_bminus-%s</sender>' % id)
    con.append('  <signal>clicked()</signal>')
    con.append('  <receiver>_expw-%s</receiver>' % id)
    con.append('  <slot>hide()</slot>')    
    con.append(' </connection>')

    return ret, con, None, None
  else:
    ret.append('<!-- /expansion -->')
    ret.append('</layout> <!-- _exp-%s-->' % id)
    return ret, None, None, None
    

header_count = {}

def generator_header(strmempos, name, id, objs, forms, wrappings, mode, args):
  if mode != "start": return None, None, None, None
  form = forms[0]
  if len(forms) == 2 and hasattr(form, "group") and form.group is not None:
    if wrappings == 2: strmempos = "!!NEW-PMEMBERPOS!!"
  if id not in header_count: header_count[id] = 0
  header_count[id] += 1
  idstr = str(id) + "-" + str(header_count[id])

  ret = []
  ret.append('<item row="%s" column="0" colspan="2">' % strmempos)
  ret.append(' <widget class="QLabel" name="_header-%s">' % idstr)
  ret.append('  <property name="font">')
  ret.append('   <font>')
  ret.append('    <italic>true</italic>')
  ret.append('   </font>')  
  ret.append('  </property>')
  ret.append('  <property name="text">')
  ret.append('   <string>%s</string>' % args[0])
  ret.append('  </property>')
  ret.append(' </widget>')
  ret.append('</item>')
  return ret, None, None, None


button_count = {}

def generator_button(strmempos,name, id, objs, forms, wrappings, mode, args):
  assert len(args) == 2, len(args)
  assert args[-1] in ("after", "before")
  if args[-1] == "before":
    if mode != "start": return None, None, None, None
  elif args[-1] == "after":
    if mode == "start": return None, None, None, None
  text = args[0]

  if id not in button_count: button_count[id] = 0
  button_count[id] += 1
  idstr = str(id) + "-" + str(button_count[id])
  
  form = forms[0]
  
  ret = []
  ret.append('<item row="%s" column="0" colspan="2">' % strmempos)
  ret.append(' <widget class="QPushButton" name="_button-%s">' % idstr)
  ret.append('  <property name="text">')
  ret.append('   <string>%s</string>' % text)
  ret.append('  </property>')
  ret.append(' </widget>')
  ret.append('</item>')

  return ret, None, None, None
 
def generator_box_layout(ori, name, id, objs, forms, wrappings, mode):

  if ori == "horizontal": ostr = "H"
  elif ori == "vertical": ostr = "V"
  else: raise ValueError(ori)
  if mode == "start": 
    ret = []
    ret.append('<layout class="Q%sBoxLayout" name="_layout-%s">' % (ostr, id))
  else:
    ret = []
    ret.append('</layout> <!-- %s -->' % id)
  return ret, None, None, None

def generator_wrap_item(strmempos, gen, *args):
  mode = args[-1]
  id = args[1]
  wrappings = args[-2]
  forms = args[-3]
  assert mode in ("start", "end")
  if wrappings == 0: #toplevel form, no outerwrap or head elements
    return gen(*args)
  ret = []
  if mode == "start" :ret.append('<item row="%s" column="0" colspan="2">' % strmempos)
  r = gen(*args)
  for rr in r[0]:
    ret.append(" " + rr)
  if mode == "end" :ret.append("</item> <!-- %s -->" % id)
  return ret, r[1], r[2], r[3]

def wrap_item(gen):
  return partial(generator_wrap_item, "!!NEW-PMEMBERPOS!!", gen)

def generator_wrap_expanded_item(gen, *args):
  mode = args[-1]
  id = args[1]
  wrappings = args[-2]
  forms = args[-3]
  assert mode in ("start", "end")
  if wrappings == 0: #toplevel form, no outerwrap or head elements
    return gen(*args)
  ret = []
  if mode == "start":
    #ret.append('<item row="!!NEW-MEMBERPOS!!" column="0" colspan="2">')
    ret.append('<item>')
    ret.append(' <widget class="QGroupBox" name="_expw-%s">' % id)
    ret.append('  <property name="__hidden_at_startup__" stdset="0">')
    ret.append('   <bool>true</bool>')
    ret.append('  </property>')      
  r = gen(*args)
  for rr in r[0]:
    ret.append("  " + rr)
  if mode == "end":
    ret.append(" </widget> <!-- _expw-%s -->" % id)
    ret.append("</item>")

  return ret, r[1], r[2], r[3]

def wrap_expanded_item(gen):
  return partial(generator_wrap_expanded_item, gen)

def generator_expander_head(name, id, objs, forms, wrappings, name2):
  ret = []
  ret.append('   <widget class="QLabel" name="_head-%s">' % id)
  ret.append('    <property name="sizePolicy">')
  ret.append('     <sizepolicy hsizetype="Expanding" vsizetype="Preferred">')
  ret.append('      <horstretch>0</horstretch>')
  ret.append('      <verstretch>0</verstretch>')
  ret.append('     </sizepolicy>')
  ret.append('    </property>')
  ret.append('    <property name="minimumSize">')
  ret.append('     <size>')
  ret.append('      <width>600</width>')
  ret.append('      <height>0</height>')
  ret.append('     </size>')
  ret.append('    </property>')
  ret.append('    <property name="font">')
  ret.append('     <font>')
  ret.append('      <weight>75</weight>')
  ret.append('      <italic>true</italic>')
  ret.append('      <bold>true</bold>')
  ret.append('     </font>')  
  ret.append('    </property>')
  ret.append('    <property name="text">')
  ret.append('     <string>%s</string>' % name2)
  ret.append('    </property>')
  ret.append('   </widget>')
  ret.append('  </item>')
  ret.append(' </layout> <!-- _exph-%s-->' % id)
  ret.append('</item>')  
  return ret, None, None, None


def generator_head(name, id, objs, forms, wrappings, name2):
  ret = []
  ret.append('<item row="!!NEW-PMEMBERPOS!!" column="0" colspan="2">')
  ret.append(' <widget class="QLabel" name="_head-%s">' % id)
  ret.append('  <property name="font">')
  ret.append('   <font>')
  ret.append('    <weight>75</weight>')
  ret.append('    <bold>true</bold>')
  ret.append('   </font>')  
  ret.append('  </property>')
  ret.append('  <property name="text">')
  ret.append('   <string>%s</string>' % name2)
  ret.append('  </property>')
  ret.append(' </widget>')
  ret.append('</item>')
  return ret, None, None, None

def generator_form_layout(name, id, objs, forms, wrappings, mode):
  if mode == "start":
    ret = []
    ret.append('<layout class="QFormLayout" name="_layout-%s">' % id)
    ret.append(' <property name="sizeConstraint">')
    ret.append('  <enum>QLayout::SetFixedSize</enum>')
    ret.append(' </property>')    
  else:
    ret = []
    ret.append('</layout> <!-- %s -->' % id)
  return ret, None, None, None
  
def generator_column_or_form_layout(name, id, objs, forms, wrappings, mode):
  form = forms[0]
  if hasattr(form, "subtype") and form.subtype == "column":
    return generator_box_layout("horizontal", name, id, objs, forms, wrappings, mode)
  else:
    return generator_form_layout(name, id, objs, forms, wrappings, mode)      

def generator_group_expander(name, id, objs, forms, wrappings, group, groupindex, mode):
  if mode in ("pre", "post"): return None
  if id is None: groupid = "group" + str(groupindex)
  else: groupid = id + "-group" + str(groupindex)
  if mode == "start":
    ret = ['<item row="!!NEW-PMEMBERPOS!!" column="0" colspan="2">']    
    ret1 = generator_expander(groupid, groupid, objs, forms, wrappings+1, mode)
    ret2 = generator_expander_head(groupid, groupid, objs, forms, wrappings+1, group)
    ret3 = wrap_expanded_item(generator_form_layout)(groupid, groupid, objs, forms, wrappings+1, mode)
    ret0 = ret + ret1[0] + ret2[0] + ret3[0]
    return ret0,ret1[1],None,None
  elif mode == "end":
    ret = ["</item>"]
    ret1 = generator_expander(groupid, groupid, objs, forms, wrappings+1, mode)    
    ret3 = wrap_expanded_item(generator_form_layout)(groupid, groupid, objs, forms, wrappings+1, mode)
    ret3[0][0] = ret3[0][0]
    return (ret3[0] + ret1[0] + ret), None,None,None
    

def init():
  from .xml import set_generator, set_element_generator, set_group_generator, set_phase_generator

  generator_header_1 = partial(generator_header, "!!NEW-MEMBERPOS!!")
  generator_button_1 = partial(generator_button, "!!NEW-MEMBERPOS!!")
  generator_header_2 = partial(generator_header, "!!NEW-MEMBERTOK!!")
  generator_button_2 = partial(generator_button, "!!NEW-MEMBERTOK!!")

  set_element_generator("textarea", generator_textarea)
  
  set_generator("outerwrap", 0, None)
  set_phase_generator("header", 0, "pre", generator_header_1)    
  set_phase_generator("button", 0, "pre", generator_button_1) 
  set_generator("group", 0, None) #never used
  set_generator("midwrap", 0, None)
  set_generator("head", 0, None) #head is already part of generated element
  set_generator("innerwrap", 0, None)

  set_generator("outerwrap", 1, None)
  set_phase_generator("header", 1, "pre", generator_header_2)    
  set_phase_generator("button", 1, "inner", generator_button_2) 
  set_generator("group", 1, None) #never used
  set_generator("midwrap", 1, None)
  set_generator("head", 1, generator_head)
  set_generator("innerwrap", 1, wrap_item(generator_form_layout))
   
  for n in range(2, 20):
    set_generator("outerwrap", n, None)
    set_phase_generator("header", n, "pre", generator_header_1)
    set_phase_generator("button", n, "inner", generator_button_1)      
    set_generator("group", n, generator_group_expander)
    set_generator("midwrap", n, wrap_item(generator_expander))
    set_generator("head", n, generator_expander_head)
    set_generator("innerwrap", n, wrap_expanded_item(generator_form_layout))

  set_generator("innerwrap", 2, wrap_expanded_item(generator_column_or_form_layout))
  
