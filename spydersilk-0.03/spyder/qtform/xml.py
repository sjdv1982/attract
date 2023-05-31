import spyder
import Spyder
from Spyder import Object, File, String, Integer, Bool, Float, Resource
from spyder.formtools.resourcewrapper import embedded_resource_filename
from functools import partial
from . import reserved

###############################################################################
# SETTINGS
###############################################################################

form_head = """<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>Form</class>"""
form_tail = """</ui>"""

form_start =""" <widget class="QScrollArea" name="Form">
  <widget class="QFrame" name="Form_contents">
   <property name="windowTitle">
    <string>Spyder.qtform</string>
   </property>"""
form_end ="""  </widget>
 </widget>"""

tabstep = 1

#if "wrappings" is between these values, the generated XML depends on "wrappings"
# therefore, if wrappings changes (due to grouping), the XML must then be regenerated
wrappings_regen_min = 1
wrappings_regen_max = 1000

###############################################################################

def indent_xml(lines):
  ret = []
  indent = 0
  for lnr, l in enumerate(lines):
    l = l.lstrip()
    pos = 0
    while 1:
      newpos = l.find("<", pos)
      if newpos != 0: break
      ch = l[newpos+1]
      if ch == "/": 
        indent -= 1
      break
    
    cindent = max(indent, 0)
    ret.append(cindent * " " + l)

    pos = 0
    while 1:
      newpos = l.find("<", pos)    
      if newpos == -1: break    
      ch = l[newpos+1]
      if ch == "/":
        if newpos == 0: pass
        else: indent -= 1
      elif ch in ("?", "!"): pass
      else: indent += 1
      pos = newpos+1  
  
  return ret
  
_generators = {}
_group_generators = {}
_ele_generators = {}
_phase_generators = {}
_tokens = ("outerwrap", "midwrap", "head", "innerwrap", "group") 
_phases = ("outer", "pre", "post", "mid", "inner")
_allowedformtypes = ("none", "hidden", "text", "password", "number", "slider", "spin", "checkbox", "option", "radio", "file")  

def set_generator(token, complexity, callback):
  """
  See the HTML generator
  """
  assert isinstance(complexity, int) and complexity >= 0
  _generators[token, complexity] = callback

def set_group_generator(grouptype, complexity, callback):
  assert isinstance(complexity, int) and complexity >= 0
  _group_generators[grouptype, complexity] = callback

def set_element_generator(element, callback):
  _ele_generators[element] = callback

def set_phase_generator(element, complexity, phase, callback):
  """
  See htmlform
  """
  assert phase in _phases, phase
  _phase_generators[element, complexity] = phase, callback
  
from . import generators
generators.init()

###########################################################################

def up2(xml, response):
  if response is None: return False
  if isinstance(response, str):
    xml.append(response)
  else:
    xml += response
  return True

def update(xml, response):
  if response is None: return False, False, False
  ret1 = up2(xml[0], response[0])
  ret2 = up2(xml[1], response[1])
  ret3 = up2(xml[2], response[2])
  ret4 = up2(xml[3], response[3])
  return ret1, ret2, ret3, ret4

def substitute_mempos(lines, exhaustive = False):
  while 1:
    ret = []
    memberpos = 0  
    replaced = False
    for l in lines:
      tokpat = "!!TOK-MEMBERPOS!!"
      tokpos = l.find(tokpat)
      if l.find("!!NEW-PMEMBERPOS!!") > -1:
        l = l.replace("!!NEW-PMEMBERPOS!!", "!!NEW-MEMBERPOS!!")      
        replaced = True
      elif l.find("!!PMEMBERPOS!!") > -1:
        l = l.replace("!!PMEMBERPOS!!", "!!MEMBERPOS!!")
        replaced = True
      elif l.find("!!NEW-MEMBERPOS!!") > -1:
        memberpos += 1
        l = l.replace("!!NEW-MEMBERPOS!!", str(memberpos-1))
        replaced = True      
      elif l.find("!!NEW-MEMBERTOK!!") > -1:
        memberpos += 1
        tokstr = "!!TOK-MEMBERPOS!!%d!!" % (memberpos - 1)
        l = l.replace("!!NEW-MEMBERTOK!!", tokstr)
        replaced = True      
      elif tokpos > -1:
        ll = l[tokpos+len(tokpat):]
        tokend = ll.find("!!")
        memberpos = int(ll[:tokend]) + 1
        tokstr = tokpat + str(memberpos-1) + "!!"
        l.index(tokstr)
        l = l.replace(tokstr, str(memberpos-1))
        replaced = True      
      elif l.find("!!MEMBERPOS!!") > -1:
        l = l.replace("!!MEMBERPOS!!", str(memberpos-1))
        replaced = True      
      ret.append(l)
    if not replaced or not exhaustive: break
    lines = ret
  return ret

def add_label(xml, n, name2):
  xml[0].append('<item row="!!NEW-MEMBERPOS!!" column="0">')
  xml[0].append('<widget class="QLabel" name="_head-%s">' % n)
  xml[0].append('<property name="text">')
  xml[0].append('<string>%s</string>' % name2) 
  xml[0].append('</property>')
  xml[0].append('</widget>')
  xml[0].append('</item>')

from .resolve import resolve
      
def gen(xml,  mode, generator, args):
  if generator is None: return
  response = generator(*args)
  if mode == "start": 
    up = update(xml, response)

  elif mode == "end": 
    up = update(xml, response)
    if response[0] is not None and not isinstance(response[0], str):
      l = response[0][0]
      dif = len(l) - len(l.lstrip())
      difpos = dif
      if l.lstrip().startswith("</"): difpos += 1
      xml[0][:] = xml[0][:-len(response[0])]
      if response[1] is not None:
        if isinstance(response[1], str): xml[1] = xml[1][:-1]
        else: xml[1][:] = xml[1][:-len(response[1])]
      if response[2] is not None:
        if isinstance(response[2], str): xml[2] = xml[2][:-1]
        else: xml[2][:] = xml[2][:-len(response[2])]
      if response[3] is not None:
        if isinstance(response[3], str): xml[3] = xml[3][:-1]
        else: xml[3][:] = xml[3][:-len(response[3])]
      update(xml, response)

def process_other_tokens(phase, xml, props, complexity, data):
  if "_othertokens" not in props: return
  name, id, objs, forms, wrappings, mode = data
  tokens = props["_othertokens"]
  if mode == "end": tokens = reversed(tokens)
  
  for token, args in tokens:
    if (token, complexity) not in _phase_generators:
      raise ValueError("Unknown token %s, complexity %s" % (token, complexity))
    tokphase, callback = _phase_generators[token, complexity]
    if tokphase != phase: continue
    cb_args = (name, id, objs, forms, wrappings, mode, args)
    gen(xml, mode, callback, cb_args)

def _write_xml(objs, forms, name=None,prename="",wrappings=0):  

  obj, form = objs[0], forms[0]

  if hasattr(form, "value"):
    defaultvalue = form.value
  elif obj is not None:
    defaultvalue = obj
  elif hasattr(form, "default"):
    defaultvalue = form.default
  else:
    defaultvalue = None
  
  complexity = 0
  mxmls = []
  mxmldata = [] #stores some values in case we have to regenerate XML
  
  n0, n = None, None
  if form.typename is None:
    typnam = "None" + "Array" * form.arraycount
  else:
    typnam = form.typename + "Array" * form.arraycount
  if name is not None:
    n0 = prename + name
    n = prename.replace(".", "-") + name

  if hasattr(form, "type"):
    if form.type is None: form.type = "none"       
    
  curr_membernames = None
  if form.arraycount > 0:
    is_array = True
    is_none = False   
    if hasattr(form, "type"):
      if form.type not in ("none","hidden"): is_array = False
      if form.type in ("none", "hidden"): is_none = True
    if is_array:
      if is_none is False:
        if not hasattr(form, "length"):
          raise TypeError("Cannot generate XML for \"%s\"(%s): member is an Array, but no length has been specified" % (n0, typnam))      
        if not isinstance(form.length, int) or form.length <= 0:
          raise TypeError("Cannot generate XML for \"%s\"(%s): member is an Array, but its length is not a positive integer" % (n0, typnam))                  
        curr_membernames = []
        curr_members = []
        curr_objs = []
        for nr in range(form.length):
          curr_membernames.append(str(nr))
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
      curr_members.append(getattr(form, membername))
      mobj = None
      if defaultvalue is not None: 
        mobj = getattr(defaultvalue, membername)
      curr_objs.append(mobj)
      
  mwrappings = wrappings + 1
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
    mcomplexities = []
    for membername, mform, mobj in zip(curr_membernames, curr_members, curr_objs):
      newprename = prename
      if name is not None: 
        newprename += name + "."
      mobjs = [mobj, defaultvalue] + objs[1:]
      mforms = [mform] + forms
      mxml, mcomplexity = _write_xml(mobjs, mforms, membername,newprename,mwrappings)
      mxmls.append(mxml)
      mxmldata.append((mobjs, mforms, newprename))
      if hasattr(mform, "group") and mform.group is not None:
        mcomplexity += 1
      mcomplexities.append(mcomplexity)
    complexity = max(mcomplexities) 
    if form.arraycount == 0 or complexity == 0: complexity +=  1  
 
  for token in _tokens:
    assert token,complexity in _generators  
  gen_group = _generators["group", complexity]
  gen_innerwrap = _generators["innerwrap", complexity]
  if name is None:
    gen_head = None
    gen_outerwrap = None
    gen_midwrap = None
  else:
    gen_head = _generators["head", complexity]    
    gen_outerwrap = _generators["outerwrap", complexity]  
    gen_midwrap = _generators["midwrap", complexity]
  xml = [],[],[],[]
  while 1:
    if name is not None:
      if form.is_default:
        if hasattr(form, "type") and form.type in ("none", "hidden"): break
      else:
        if defaultvalue is None and hasattr(form, "type") and form.type in ("none", "hidden"):
          raise TypeError("Cannot generate XML for \"%s\"(%s): member is not optional, but type is specified as \"%s\"" % (n0, typnam, form.type))

    tokdata = (name, n, objs, forms, wrappings, "start")

    process_other_tokens("outer", xml, form._props, complexity, tokdata)

    gen_args = (name, n, objs, forms, wrappings, "start")
    gen(xml, "start", gen_outerwrap, gen_args)

    process_other_tokens("pre", xml, form._props, complexity, tokdata)

    gen(xml, "start", gen_midwrap, gen_args)

    process_other_tokens("mid", xml, form._props, complexity, tokdata)
    
    xml0 = list(xml[0]), list(xml[1]), list(xml[2]), list(xml[3])
    xml = [], [], [], []

    name2 = name
    if name2 is not None and name2.endswith("_"):
      for r in reserved:
        if name2 == r + "_": 
          name2 = r
          break
    if hasattr(form, "name"): name2 = form.name
    if gen_head is not None: 
      val = gen_head(name, n, objs, forms, wrappings, name2)
      update(xml, val)

    process_other_tokens("post", xml, form._props, complexity, tokdata)

    gen(xml, "start", gen_innerwrap, gen_args)

    process_other_tokens("inner", xml, form._props, complexity, tokdata)
    
    curr_wrappings = mwrappings
    has_widget_type = False
    if hasattr(form, "type") and form.type not in (None, "none", "hidden"): has_widget_type = True
    if complexity > 0 and not has_widget_type:
      done_membernames = set()
      groupindex = 0
      resolve(
       xml,
       curr_wrappings, 
       done_membernames,
       mxmls,
       mxmldata,
       mwrappings,
       curr_membernames,
       curr_members,
       mcomplexities,
       groupindex,
       name,
       n,
       objs,
       forms,
       wrappings,
       complexity,
       gen_group,
       toplevel = True
      )
    else: #complexity 0, or form.type
      try:
        typ = getattr(Spyder, form.typename)
      except AttributeError:
        typ = Spyder.String
      formtype = "text"
      if issubclass(typ, File):
        if form.arraycount > 0:
          raise TypeError("Cannot generate XML for \"%s\"(%s): array of files" % (n0, typnam))
        if not hasattr(form, "file"):
          raise TypeError("Cannot generate XML for \"%s\"(%s): a file type (file attribute) must be specified" % (n0, typnam))
        if hasattr(form, "options"):
          raise TypeError("Cannot generate XML for \"%s\"(%s): file elements cannot have options" % (n0, typnam))
        formtype = "file"
      elif issubclass(typ, Resource):
        formtype = "file"
      elif issubclass(typ, Integer) or issubclass(typ, Float):
        formtype = "spin"
        
      elif issubclass(typ, Bool):
        formtype = "checkbox"
       
      if hasattr(form, "options"):
        assert formtype != "file"
        formtype = "option"
        options = form.options      
        optiontitles = options
        if hasattr(form, "optiontitles"):
          optiontitles = form.optiontitles
          assert len(options) == len(optiontitles), (n, len(options), len(optiontitles))
     
      if hasattr(form, "type"):
        if form.type not in ("none","hidden"):
          oldformtype = formtype
          formtype = form.type
          if hasattr(form._typetree, "typemapped") and form._typetree.typemapped == True:
            if formtype not in _allowedformtypes and formtype not in _ele_generators:
              formtype = oldformtype
        
      assert formtype in _allowedformtypes or formtype in _ele_generators, formtype

      formsubtype = None
      if hasattr(form, "subtype"):
        formsubtype = form.subtype
      if formsubtype is None:
        if issubclass(typ, Integer):
          formsubtype = "int"
        elif issubclass(typ, Float):
          formsubtype = "float"

      expand = True
      if hasattr(form, "expand"): expand = form.expand
                  
      if formtype in _ele_generators:
        gener = _ele_generators[formtype]
        ret = gener(name2, n, objs, forms, wrappings, defaultvalue) 
        update(xml, ret)
      else: #standard types
        if formtype == "option":
          assert hasattr(form, "options"), n
          
          add_label(xml, n, name2)
          xml[0].append('<item row="!!MEMBERPOS!!" column="1">')
          xml[0].append('<widget class="QComboBox" name="_widget-%s">' % n)
          if defaultvalue is not None:
            assert defaultvalue in options, (n, defaultvalue, options)
            currentIndex = options.index(defaultvalue)
          else:
            currentIndex = -1
          xml[0].append('<property name="currentIndex">')
          xml[0].append('<number>%d</number>' % currentIndex)
          xml[0].append('</property>')
          
            
                      
          for onr in range(len(options)):
            opt, optitle = options[onr], optiontitles[onr]           
            xml[0].append('<item>')
            xml[0].append('<property name="text">')
            xml[0].append('<string>%s</string>' % optitle)
            xml[0].append('</property>')
            xml[0].append('</item>')         

          xml[0].append('</widget>')
          xml[0].append('</item>')

        elif formtype == "radio":
          assert hasattr(form, "options"), n
          if defaultvalue is not None:
            assert defaultvalue in options, (n, defaultvalue, options)          

          xml[3].append('<buttongroup name="_widget-%s"/>' % n)
          add_label(xml, n, name2)
          xml[0].append('<item row="!!MEMBERPOS!!" column="1">')
          xml[0].append('<layout class="QVBoxLayout" name="_box-%s">' % n)

          for onr in range(len(options)):
            opt, optitle = options[onr], optiontitles[onr]
            oid = "_widget-%s-%d" % (n, onr)
            if onr == 0: oid = "_widget-%s" % n
            
            xml[0].append('<item>')
            xml[0].append('<widget class="QRadioButton" name="%s">' % oid)
            xml[0].append('<property name="text">')
            xml[0].append('<string>%s</string>' % optitle)
            xml[0].append('</property>')
            if defaultvalue is not None and opt == defaultvalue:
              xml[0].append('<property name="checked">')
              xml[0].append('<bool>true</bool>')
              xml[0].append('</property>')
            xml[0].append('<property name="buttonGroup">')
            xml[0].append('<string notr="true">"_widget-%s"</string>' % n)
            xml[0].append('</property>')
            xml[0].append('</widget>')
            xml[0].append('</item>')
            
          xml[0].append('</layout>')
          xml[0].append('</item>')

        elif formtype == "file":
          defaultstr = ""
          if defaultvalue is not None:
            if isinstance(defaultvalue, File):
              defaultstr = defaultvalue.name
            elif isinstance(defaultvalue, Resource):
              if defaultvalue.filename is not None:
                defaultstr = defaultvalue.filename
              else:
                defaultstr = embedded_resource_filename
            elif isinstance(defaultvalue, str):
              defaultstr = defaultvalue

          add_label(xml, n, name2)
          xml[0].append('<item row="!!MEMBERPOS!!" column="1">')
          xml[0].append('<layout class="QHBoxLayout" name="_box-%s">' % n)
          xml[0].append('<item>')
          xml[0].append('<widget class="QLineEdit" name="_widget-%s">' % n)

          xml[0].append('<property name="readOnly">')
          xml[0].append('<bool>true</bool>') 
          xml[0].append('</property>')
          if defaultstr != "":
            xml[0].append('<property name="text">')
            xml[0].append('<string>%s</string>' % defaultstr) 
            xml[0].append('</property>')
          xml[0].append('</widget>')          
          xml[0].append('</item>')
          
          xml[0].append('<item>')
          xml[0].append('<widget class="QPushButton" name="_selector-%s">' % n)

          xml[0].append('<property name="text">')
          xml[0].append('<string>Browse...</string>')
          xml[0].append('</property>')

          xml[0].append('<property name="__file_select__" stdset="0">')
          xml[0].append('<bool>true</bool>') 
          xml[0].append('</property>')

          if hasattr(form, "folder") and form.folder == True:
            xml[0].append('<property name="__file_select_folder__" stdset="0">')
            xml[0].append('<bool>true</bool>') 
            xml[0].append('</property>')

          xml[0].append('</widget>')          
          xml[0].append('</item>')
            
          xml[0].append('</layout>')
          xml[0].append('</item>')

        elif formtype == "checkbox":

          xml[0].append('<item row="!!NEW-MEMBERPOS!!" column="0" colspan="2">')
          xml[0].append('<widget class="QCheckBox" name="_widget-%s">' % n)
          xml[0].append('<property name="text">')
          xml[0].append('<string>%s</string>' % name2)
          xml[0].append('</property>')
          if defaultvalue is not None:
            xml[0].append('<property name="checked">')
            val = str(bool(defaultvalue)).lower()
            xml[0].append('<bool>%s</bool>' % val)
            xml[0].append('</property>')
          
          xml[0].append('</widget>')
          xml[0].append('</item>')

        elif formtype == "spin":
          step = 1
          if hasattr(form, "step"): step = form.step
          elif hasattr(form, "range"): 
            step = form.range/10.0
          step = min(step, 1)  
          if formsubtype == "int": step = int(step+0.999)

          digits = None
          if formsubtype == "float":
            if hasattr(form, "digits"): 
              digits = form.digits          
              digits = max(digits, 2)  
              step = 10**-digits
            elif step < 1 and step > 0:
              digits = 1
              while 10**-digits > step: digits += 1
              digits = max(digits, 2)  
            elif defaultvalue is not None:
              trail = str(float(defaultvalue)-int(defaultvalue)).rstrip("0")
              if trail.endswith("."): 
                digits = 0
              else:
                digits = len(trail) - 2 #len("0.")  
              digits = max(digits, 2)    
              step = 10**-digits
            else:
              digits = 2
              step = 0.01  

          minimum = -999999
          if hasattr(form, "min"): minimum = form.min
          maximum = 999999
          if hasattr(form, "max"): maximum = form.max

          if formsubtype == "float": 
            w = "QDoubleSpinBox"
            printprop = '<double>%f</double>'
          else: 
            w = "QSpinBox"
            printprop = '<number>%d</number>'

          add_label(xml, n, name2)
          xml[0].append('<item row="!!MEMBERPOS!!" column="1">')
          xml[0].append('<widget class="%s" name="_widget-%s">' % (w,n))

          if digits is not None:
            xml[0].append('<property name="decimals">')
            xml[0].append('<number>%d</number>' % digits)
            xml[0].append('</property>')
          
          xml[0].append('<property name="minimum">')
          xml[0].append(printprop % minimum)
          xml[0].append('</property>')

          xml[0].append('<property name="maximum">')
          xml[0].append(printprop % maximum) 
          xml[0].append('</property>')

          xml[0].append('<property name="singleStep">')
          xml[0].append(printprop % step)
          xml[0].append('</property>')

          if defaultvalue is not None:
            xml[0].append('<property name="value">')
            xml[0].append(printprop % defaultvalue)
            xml[0].append('</property>')

          xml[0].append('</widget>')
          xml[0].append('</item>')

        elif formtype == "text" or formtype == "password":
          add_label(xml,n, name2)
          defaultstr = ""
          if defaultvalue is not None:
            defaultstr = str(defaultvalue)

          xml[0].append('<item row="!!MEMBERPOS!!" column="1">')
          xml[0].append('<widget class="QLineEdit" name="_widget-%s">' % n)
          if defaultstr != "":
            xml[0].append('<property name="text">')
            xml[0].append('<string>%s</string>' % defaultstr) 
            xml[0].append('</property>')
          xml[0].append('</widget>')
          xml[0].append('</item>')

        else:        
          raise Exception(formtype)

      ### ENDIF standard types    
    ### ENDIF complexity == 0
      
    xml = xml0[0] + xml[0], xml0[1] + xml[1], xml0[2] + xml[2], xml0[3] + xml[3]
    tokdata = (name, n, objs, forms, wrappings, "end")
    process_other_tokens("inner",  xml, form._props, complexity, tokdata)
    
    gen_args = (name, n, objs, forms, wrappings, "end")
    
    gen(xml, "end", gen_innerwrap, gen_args)

    process_other_tokens("post", xml, form._props, complexity, tokdata)
    
    gen(xml, "end", gen_midwrap, gen_args)

    process_other_tokens("mid", xml, form._props, complexity, tokdata)
    process_other_tokens("pre", xml, form._props, complexity, tokdata)
    
    gen(xml, "end", gen_outerwrap, gen_args)

    process_other_tokens("outer", xml, form._props, complexity, tokdata)
    
    if complexity > 0:
      newxml = substitute_mempos(xml[0])
      xml = (newxml,) + xml[1:]
      
    break

  return xml, complexity

from . import generators

def xml(obj=None, form=None, head = None, tail = None, start = None, end = None, indent1 = 3, indent2 = 2):
  """
  obj can be: A Spyder class, a Spyder object, or None
    if obj is None, then form must not be None
  form is a spyderform object
    if form is None, it is extracted from the Spyder obj
  
  xml[0]: main widget elements
  xml[1]: connections
  xml[2]: resources
  xml[3]: buttongroups
  
  """
  generators.button_count.clear()
  if head is None: head = form_head
  if tail is None: tail = form_tail
  if start is None: start = form_start
  if end is None: end = form_end
  
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
  from ..formtools import arraymanager_dynamic
  arraymanager_dynamic.add_buttons(form)
   
  xml, complexity = _write_xml([obj], [form])
  
  xml = [indent_xml(x) for x in xml]
  
  ret = head + "\n"
  ret += start + "\n"
  
  newxml = substitute_mempos(xml[0], True)  
  for l in newxml:
    ret += indent1 * " " + l + "\n"
  ret += end + "\n"
  ret += (indent2-1) * " " + "<connections>\n"  
  for l in xml[1]:
    ret += indent2 * " " + l + "\n"
  ret += (indent2-1) * " " + "</connections>\n"  
  ret += (indent2-1) * " " + "<resources>\n"  
  for l in xml[2]:
    ret += indent2 * " " + l + "\n"
  ret += (indent2-1) * " " + "</resources>\n"  
  ret += (indent2-1) * " " + "<buttongroups>\n"  
  for l in xml[3]:
    ret += indent2 * " " + l + "\n"
  ret += (indent2-1) * " " + "</buttongroups>\n"  
  ret += tail + "\n"
    
  return ret
    
  

