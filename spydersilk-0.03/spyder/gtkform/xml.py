import spyder
import Spyder
from Spyder import Object, File, String, Integer, Bool, Float, Resource
from spyder.formtools.resourcewrapper import embedded_resource_filename
from functools import partial

###############################################################################
# SETTINGS
###############################################################################

#TODO: groups, check if they work

#TODO: document .folder, .subtype="column"
#TODO: styles

form_head = """<?xml version="1.0" encoding="UTF-8"?>
<interface>
  <!-- interface-requires gtk+ 3.0 -->"""
form_tail = "</interface>"

form_start ="""  <object class="GtkWindow" id="window1">
    <property name="can_focus">False</property>
    <property name="has_resize_grip">False</property>
    <child>
      <object class="GtkScrolledWindow" id="scrolledwindow1">
        <property name="visible">True</property>
        <property name="can_focus">True</property>
        <property name="shadow_type">in</property>
        <property name="width_request">800</property>
        <property name="height_request">600</property>
        <child>
          <object class="GtkViewport" id="viewport1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <child>"""   
    
form_end ="""        
            </child>
          </object>        
        </child>
      </object>
    </child>      
  </object>"""

tabstep = 2

#if "wrappings" is between these values, the generated XML depends on "wrappings"
# therefore, if wrappings changes (due to grouping), the XML must then be regenerated
wrappings_regen_min = 1
wrappings_regen_max = 1000

###############################################################################

_generators = {}
_group_generators = {}
_ele_generators = {}
_phase_generators = {}
_tokens = ("outerwrap", "midwrap", "head", "innerwrap", "group") 
_phases = ("outer", "pre", "post", "mid", "inner")
_allowedformtypes = ("none", "text", "password", "number", "slider", "spin", "checkbox", "option", "radio", "file")  

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

def up2(xml, response, tabstr):
  if response is None: return False
  if isinstance(response, str):
    xml.append(tabstr + response)
  else:
    for l in response:
      xml.append(tabstr + l)
  return True

def update(xml, response, t, dif=0):
  if response is None: return False, False
  ret1 = up2(xml[0], response[0], (t.currtabs + dif) * " ")
  ret2 = up2(xml[1], response[1], (t.currtabs2 + dif) * " ")
  return ret1, ret2

from .resolve import resolve
      
def gen(xml, t, mode, generator, args):
  if generator is None: return
  response = generator(*args)
  if mode == "start": 
    up = update(xml, response, t, 0)
    if up[0]: t.incr()
    if up[1]: t.incr2()      

    if response[0] is not None and not isinstance(response[0], str):
      l = response[0][-1]
      dif = len(l) - len(l.lstrip())
      difpos = dif / 2
      if l.find("</") > -1: difpos -= 1
      if difpos > 0: 
        for n in range(int(difpos)): t.incr()
  elif mode == "end": 
    up = update(xml, response, t, -tabstep)
    if up[0]: t.decr()
    if up[1]: t.decr2()      
    if response[0] is not None and not isinstance(response[0], str):
      l = response[0][0]
      dif = len(l) - len(l.lstrip())
      difpos = dif / 2
      if l.lstrip().startswith("</"): difpos += 1
      if difpos > 0:         
        xml[0][:] = xml[0][:-len(response[0])]
        if response[1] is not None:
          if isinstance(response[1], str): xml[1] = xml[1][:-1]
          else: xml[1][:] = xml[1][:-len(response[1])]
        for n in range(int(difpos)): t.decr()        
        update(xml, response, t)

def process_other_tokens(phase, t, xml, props, complexity, data):
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
    gen(xml, t, mode, callback, cb_args)

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

  class tabclass:
    def __init__(self):
     self.currtabs = 0
     self.tabstr = ""    
     self.currtabs2 = 0
     self.tabstr2 = ""    
    def incr(self):
      self.currtabs += tabstep
      self.tabstr = self.currtabs * " "
    def decr(self):
      self.currtabs -= tabstep
      self.tabstr = self.currtabs * " "
    def incr2(self):
      self.currtabs2 += tabstep
      self.tabstr2 = self.currtabs2 * " "
    def decr2(self):
      self.currtabs2 -= tabstep
      self.tabstr2 = self.currtabs2 * " "
  t = tabclass()
  
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
      if form.type not in ("none",): is_array = False
      if form.type == "none": is_none = True
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
        assert m in curr_membernames, m        
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
  xml = [],[]
  while 1:
    if name is not None:
      if form.is_default:
        if hasattr(form, "type") and form.type == "none": break
      else:
        if hasattr(form, "type") and form.type == "none":
          raise TypeError("Cannot generate XML for \"%s\"(%s): member is not optional, but type is specified as \"none\"" % (n0, typnam))

    tokdata = (name, n, objs, forms, wrappings, "start")

    process_other_tokens("outer", t, xml, form._props, complexity, tokdata)

    gen_args = (name, n, objs, forms, wrappings, "start")
    gen(xml, t, "start", gen_outerwrap, gen_args)

    process_other_tokens("pre", t, xml, form._props, complexity, tokdata)

    gen(xml, t, "start", gen_midwrap, gen_args)

    process_other_tokens("mid", t, xml, form._props, complexity, tokdata)
    
    xml0 = list(xml[0]), list(xml[1])
    xml = [], []

    name2 = name
    if hasattr(form, "name"): name2 = form.name
    if gen_head is not None: 
      val = gen_head(name, n, objs, forms, wrappings, name2)
      update(xml, val, t)

    process_other_tokens("post", t, xml, form._props, complexity, tokdata)

    gen(xml, t, "start", gen_innerwrap, gen_args)

    process_other_tokens("inner", t, xml, form._props, complexity, tokdata)
    
    curr_wrappings = mwrappings
    has_widget_type = False
    if hasattr(form, "type") and form.type not in (None, "none"): has_widget_type = True
    if complexity > 0 and not has_widget_type:
      done_membernames = set()
      groupindex = 0
      resolve(
       xml,
       t,
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
      typ = getattr(Spyder, form.typename)
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
        if form.type not in ("none",):
          formtype = form.type
        
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
        ret = gener(name, n, objs, forms, wrappings, defaultvalue) 
        update(xml, ret, t)
      else: #standard types

        xml[0].append(t.tabstr + '<object class="GtkLabel" id="_head-%s">' % n)
        t.incr()
        xml[0].append(t.tabstr + '<property name="visible">True</property>')
        xml[0].append(t.tabstr + '<property name="hexpand">%s</property>' % expand)
        xml[0].append(t.tabstr + '<property name="vexpand">%s</property>' % expand)
        xml[0].append(t.tabstr + '<property name="can_focus">False</property>')        
        #xml[0].append(t.tabstr + '<property name="label" translatable="yes">&lt;span foreground="blue" size="medium"&gt;%s&lt;/span&gt;</property>' % name2) #TODO
        xml[0].append(t.tabstr + '<property name="label" translatable="yes">%s</property>' % name2) 
        #xml[0].append(t.tabstr + '<property name="use_markup">True</property>') #TODO
        t.decr()
        xml[0].append(t.tabstr + '</object>')        
        xml[0].append(t.tabstr + '<packing>')
        t.incr()
        xml[0].append(t.tabstr + '<property name="expand">%s</property>' % expand)        
        xml[0].append(t.tabstr + '<property name="fill">%s</property>' % expand)        
        xml[0].append(t.tabstr + '<property name="left_attach">0</property>')
        xml[0].append(t.tabstr + '<property name="top_attach">!!NEW-MEMBERPOS!!</property>')
        xml[0].append(t.tabstr + '<property name="width">1</property>')
        xml[0].append(t.tabstr + '<property name="height">1</property>')        
        t.decr()
        xml[0].append(t.tabstr + '</packing>')
        t.decr()
        xml[0].append(t.tabstr + '</child>')        
        xml[0].append(t.tabstr + '<child>')                
        t.incr()
        if formtype == "option":
          assert hasattr(form, "options"), n
          if defaultvalue is not None:
            assert defaultvalue in options, (n, defaultvalue, options)          
          t.incr()
          xml[0].append(t.tabstr + '<object class="GtkComboBoxText" id="_widget-%s">' % n)
          xml[0].append(t.tabstr + '<property name="visible">True</property>')          
          xml[0].append(t.tabstr + '<property name="can_focus">True</property>')
          xml[0].append(t.tabstr + '<property name="hexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="vexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<items>')
          t.incr()
          for opt, optitle in zip(options, optiontitles):
            #TODO: retrieval mechanism: match option <=> optiontitle
            it = '<item translatable="yes">'
            op = str(optitle)
            xml[0].append(t.tabstr + it + op + "</item>")
          t.decr()
          xml[0].append(t.tabstr + '</items>')
          if defaultvalue is not None:
            for onr,opt in enumerate(options):
              #Does not currently work (bug in Glade): 
              # http://stackoverflow.com/questions/14912210/set-gtk-comboboxtext-default-item 
              if opt == defaultvalue:
                xml[0].append(t.tabstr + '<property name="active">%s</property>' % onr)          
                break
          t.decr()
          xml[0].append(t.tabstr + '</object>')        
        elif formtype == "radio":
          assert hasattr(form, "options"), n
          if defaultvalue is not None:
            assert defaultvalue in options, (n, defaultvalue, options)          
          xml[0].append(t.tabstr + '<object class="GtkBox" id="box_%s">' % n)
          t.incr()
          xml[0].append(t.tabstr + '<property name="visible">True</property>')
          xml[0].append(t.tabstr + '<property name="can_focus">False</property>')
          xml[0].append(t.tabstr + '<property name="orientation">vertical</property>')          
          for onr in range(len(options)):
            opt, optitle = options[onr], optiontitles[onr]
            oid = "_widget-%s-%d" % (n, onr)
            if onr == 0: oid = "_widget-%s" % n
            xml[0].append(t.tabstr + '<child>')
            t.incr()            
            xml[0].append(t.tabstr + '<object class="GtkRadioButton" id="%s">' % oid)
            t.incr()
            xml[0].append(t.tabstr + '<property name="label" translatable="yes">%s</property>' % optitle)
            xml[0].append(t.tabstr + '<property name="visible">True</property>')
            xml[0].append(t.tabstr + '<property name="can_focus">False</property>')
            xml[0].append(t.tabstr + '<property name="receives_default">False</property>')
            xml[0].append(t.tabstr + '<property name="use_action_appearance">False</property>')
            xml[0].append(t.tabstr + '<property name="xalign">0</property>')
            xml[0].append(t.tabstr + '<property name="draw_indicator">True</property>')
            if onr > 0:
              xml[0].append(t.tabstr + '<property name="group">%s</property>' % ("_widget-%s" % n))
            if defaultvalue is not None and opt == defaultvalue:
              xml[0].append(t.tabstr + '<property name="active">True</property>')
            t.decr()
            xml[0].append(t.tabstr + '</object>')
            xml[0].append(t.tabstr + '<packing>')
            t.incr()
            xml[0].append(t.tabstr + '<property name="expand">%s</property>' % expand)        
            xml[0].append(t.tabstr + '<property name="fill">%s</property>' % expand)        
            xml[0].append(t.tabstr + '<property name="position">%d</property>' % onr)
            t.decr()
            xml[0].append(t.tabstr + '</packing>')            
            t.decr()
            xml[0].append(t.tabstr + '</child>')
          t.decr()
          xml[0].append(t.tabstr + '</object>')        
        elif formtype == "file":
          xml[0].append(t.tabstr + '<object class="GtkFileChooserButton" id="_widget-%s">' % n)
          t.incr()
          xml[0].append(t.tabstr + '<property name="visible">True</property>')          
          xml[0].append(t.tabstr + '<property name="can_focus">True</property>')
          xml[0].append(t.tabstr + '<property name="hexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="vexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="orientation">vertical</property>')
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
            xml[0].append(t.tabstr + '<property name="active">%s</property>' % defaultstr)
          if hasattr(form, "folder") and form.folder == True:
            xml[0].append(t.tabstr + '<property name="action">select-folder</property>')
          t.decr()
          xml[0].append(t.tabstr + '</object>')
        elif formtype == "checkbox":
          xml[0].append(t.tabstr + '<object class="GtkCheckButton" id="_widget-%s">' % n)
          t.incr()
          xml[0].append(t.tabstr + '<property name="visible">True</property>')          
          xml[0].append(t.tabstr + '<property name="can_focus">True</property>')
          xml[0].append(t.tabstr + '<property name="hexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="vexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="receives_default">True</property>')
          xml[0].append(t.tabstr + '<property name="use_action_appearance">True</property>')
          xml[0].append(t.tabstr + '<property name="use_stock">True</property>')
          if defaultvalue is not None:
           xml[0].append(t.tabstr + '<property name="active">%s</property>' % str(bool(defaultvalue)))
          xml[0].append(t.tabstr + '<property name="draw_indicator">True</property>')           
          t.decr()
          xml[0].append(t.tabstr + '</object>')
        elif formtype == "spin":
          step = 1
          if hasattr(form, "step"): step = form.step
          elif hasattr(form, "range"): 
            step = form.range/10.0
          if formsubtype == "int": step = int(step+0.999)

          digits = 0
          if hasattr(form, "digits"): digits = form.digits
          elif formsubtype == "float":
            digits = 1
            if step < 1 and step > 0:
              while 10**-digits > step: digits += 1

          #spinbutton
          ####
          xml[0].append(t.tabstr + '<object class="GtkSpinButton" id="_widget-%s">' % n)
          t.incr()
          xml[0].append(t.tabstr + '<property name="visible">True</property>')
          xml[0].append(t.tabstr + '<property name="hexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="vexpand">%s</property>' % expand)
          
          if hasattr(form, "length"):
            xml[0].append(t.tabstr + '<property name="max_length">%d</property>' % form.length)
            xml[0].append(t.tabstr + '<property name="width_chars">%d</property>' % form.length)
          xml[0].append(t.tabstr + '<property name="can_focus">True</property>')
          xml[0].append(t.tabstr + '<property name="adjustment">%s</property>' % n)
          xml[0].append(t.tabstr + '<property name="snap_to_ticks">False</property>')
          xml[0].append(t.tabstr + '<property name="numeric">True</property>')            
          xml[0].append(t.tabstr + '<property name="digits">%d</property>' % digits)          
          t.decr()
          xml[0].append(t.tabstr + '</object>')
          ####
         
          #adjustment
          ####
          xml[1].append(t.tabstr2 + '<object class="GtkAdjustment" id="%s">' % n)
          t.incr2()                                
          minimum = -999999
          if hasattr(form, "min"): minimum = form.min
          xml[1].append(t.tabstr2 + '<property name="lower">%s</property>' % minimum)
          maximum = 999999
          if hasattr(form, "max"): maximum = form.max
          xml[1].append(t.tabstr2 + '<property name="upper">%s</property>' % maximum)
          if defaultvalue is not None:
            xml[1].append(t.tabstr2 + '<property name="value">%s</property>' % defaultvalue)
          
          xml[1].append(t.tabstr2 + '<property name="step_increment">%s</property>' % step)
                        
              
          t.decr2()
          xml[1].append(t.tabstr2 + '</object>')
          ####

        elif formtype == "text" or formtype == "password":
          defaultstr = ""
          if defaultvalue is not None:
            defaultstr = str(defaultvalue)
          xml[0].append(t.tabstr + '<object class="GtkEntry" id="_widget-%s">' % n)
          t.incr()
          xml[0].append(t.tabstr + '<property name="visible">True</property>')          
          xml[0].append(t.tabstr + '<property name="can_focus">True</property>')
          xml[0].append(t.tabstr + '<property name="hexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="vexpand">%s</property>' % expand)
          xml[0].append(t.tabstr + '<property name="invisible_char">\xe2\x80\xa2</property>')
          if defaultvalue is not None:
            xml[0].append(t.tabstr + '<property name="text">%s</property>' % defaultvalue)
          if hasattr(form, "length"):
            xml[0].append(t.tabstr + '<property name="max_length">%d</property>' % form.length)
            xml[0].append(t.tabstr + '<property name="width_chars">%d</property>' % form.length)
          if formtype == "password":
            xml[0].append(t.tabstr + '<property name="visibility">False</property>')
          t.decr()
          xml[0].append(t.tabstr + '</object>')
        else:        
          raise Exception(formtype)


        xml[0].append(t.tabstr + '<packing>')
        t.incr()
        xml[0].append(t.tabstr + '<property name="expand">%s</property>' % expand)        
        xml[0].append(t.tabstr + '<property name="fill">%s</property>' % expand)        
        xml[0].append(t.tabstr + '<property name="left_attach">1</property>')
        xml[0].append(t.tabstr + '<property name="top_attach">!!MEMBERPOS!!</property>')
        xml[0].append(t.tabstr + '<property name="width">1</property>')
        xml[0].append(t.tabstr + '<property name="height">1</property>')        
        t.decr()
        xml[0].append(t.tabstr + '</packing>')
        t.decr()
      ### ENDIF standard types    
    ### ENDIF complexity == 0
      
    xml = xml0[0] + xml[0], xml0[1] + xml[1]
    tokdata = (name, n, objs, forms, wrappings, "end")
    process_other_tokens("inner", t, xml, form._props, complexity, tokdata)
    
    gen_args = (name, n, objs, forms, wrappings, "end")
    
    gen(xml, t, "end", gen_innerwrap, gen_args)

    process_other_tokens("post", t, xml, form._props, complexity, tokdata)

    newxml = []
    memberpos = 0
    if complexity > 0:
      for l in xml[0]:
        if l.find("!!NEW-MEMBERPOS!!") > -1:
          memberpos += 1
          l = l.replace("!!NEW-MEMBERPOS!!", str(memberpos-1))
        elif l.find("!!MEMBERPOS!!") > -1:
          l = l.replace("!!MEMBERPOS!!", str(memberpos-1))
        newxml.append(l)
      xml = (newxml, xml[1])
    
    gen(xml, t, "end", gen_midwrap, gen_args)

    process_other_tokens("mid", t, xml, form._props, complexity, tokdata)
    process_other_tokens("pre", t, xml, form._props, complexity, tokdata)
    
    gen(xml, t, "end", gen_outerwrap, gen_args)

    process_other_tokens("outer", t, xml, form._props, complexity, tokdata)
    
      
    break

  return xml, complexity


def xml(obj=None, form=None, head = None, tail = None, start = None, end = None, indent1 = 14, indent2 = 2):
  """
  obj can be: A Spyder class, a Spyder object, or None
    if obj is None, then form must not be None
  form is a spyderform object
    if form is None, it is extracted from the Spyder obj
  """
  
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
   
  xml, complexity = _write_xml([obj], [form])
  
  ret = head + "\n"
  ret += start + "\n"
  for l in xml[0]:
    l = l.replace("!!NEW-MEMBERPOS!!", "0")
    l = l.replace("!!MEMBERPOS!!", "0")
    ret += indent1 * " " + l + "\n"
  ret += end + "\n"
  for l in xml[1]:
    ret += indent2 * " " + l + "\n"
  ret += tail + "\n"
  return ret
    
  

