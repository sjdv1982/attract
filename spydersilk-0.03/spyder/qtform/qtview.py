from functools import partial
import weakref 
import sys
from .. import validvar2
python3 = (sys.version_info[0] == 3)

class pseudometaobject(object):
  @staticmethod
  def className(): return "pseudowidget"
  
class pseudowidget(object):
  def __init__(self, widgets):
    self._widgets = widgets
  @staticmethod
  def metaObject(): return pseudometaobject
  def show(self):
    for w in self._widgets: w.show()
  def hide(self):
    for w in self._widgets: w.hide()

def collect_widgets(widgets):
  ret = []
  for w in widgets:
    if isinstance(w, tuple): 
      wl = collect_widgets(w)
      ret += wl
    else:  
      classname = w.metaObject().className()
      if classname == "QWidget": continue
      ret.append(w)
  return ret

def find_common_ancestor(widgets):
  parentchains = []
  for w in widgets:
    parentchain = []
    ww = w
    while ww is not None:
      parentchain.append(ww)
      ww = ww.parent()
    parentchain.reverse()  
    parentchains.append(parentchain)
  common = None
  for count in range(len(parentchains[0])):
    ok = True
    w = parentchains[0][count]
    for p in parentchains[1:]:
      if len(p) <= count or p[count] is not w:
        ok = False
        break
    if not ok: break
    common = w
  if common is None:
    raise ValueError("Widgets have no common ancestor")
  return common    

def make_pseudowidget(wrapwidgets, ancestor):
  widgets = []
  widgetids = set()
  for w in wrapwidgets:
    ww = w
    while ww.parent() is not ancestor: ww = ww.parent()
    if ww in widgetids: continue
    widgets.append(ww)
    widgetids.add(id(ww))
  return pseudowidget(widgets)  
  
def if_validval(func, typ, value):
  if value == None: return
  val = typ(value)
  return func(val)

def validval(func, typ):
  return partial(if_validval, func, typ)

def getpropfunc(getprop, filetype):
  return (getprop(), filetype)
  
def lineeditval(widget, value):
  assert value is not None
  val = str(value)
  curr = widget.text()
  if curr != val:
    pos = widget.cursorPosition()
    widget.setText(val)
    widget.setCursorPosition(pos)

def texteditval(widget, value):
  assert value is not None
  val = str(value)
  curr = widget.toPlainText()
  if curr != val:
    cursor = widget.textCursor()
    widget.setPlainText(val)
    widget.setTextCursor(cursor)

def selector_eval(func, func_args, target):
  result = func(*func_args)
  target(result)

def de_unicode(func):
  val = func()
  return str(val)
  
def string_getter(func):
  if python3: return func
  return partial(de_unicode, func)

def wrap_resource(getprop, setprop, form, formtype):
  from spyder.formtools.resourcewrapper import ResourceWrapper
  typename = form.typename + form.arraycount * "Array"
  typename2 = "Resource"+typename
  getprop2 = getprop
  if formtype == "file": 
    getprop2 = partial(getpropfunc, getprop, typename)
  w = ResourceWrapper(getprop2, setprop, typename2, formtype)
  return w.getprop, w.setprop
  
def openfilename(func):
  from . import anyQt
  if anyQt._qt == "PySide":
    return func()[0]
  elif anyQt._qt == "PyQt4":
    return func()    
  else: raise Exception
  
def build_widgets(qtui):
  from .anyQt.QtCore import QObject
  ret = {}
  widgets = qtui.findChildren(QObject)
  for w in widgets:
    hid = w.property("__hidden_at_startup__")
    if hid: w.hide()  
    try:
      name = str(w.objectName())
    except TypeError:
      continue
    if len(name) == 0: continue
    if name.startswith("_"):
      try:
        pos = name.index("-")
        n1, n2 = name[:pos], name[pos+1:]
        ret[n1,n2] = w
      except ValueError:
        raise ValueError("Malformed widget name %s" % name)
    else:
      ret[None,name] = w
  return ret

def trim_widgets(widgets, trim):
  ret = {}
  for cat,name in widgets:
    if name.startswith(trim):
      name2 = name[len(trim):]
      ret[cat,name2] = widgets[cat,name]
  return ret

class buttonwrapper(object):
  def __init__(self, widget):
    self.widget = widget
    self._listeners = []
    self._listen_handler = None
    self._blockedcount = 0
  def _listener(self, *args):
    #print("LISTEN!", self, self._blockedcount, self._value_callback())
    if self._blockedcount: return
    for callback in self._listeners:
      callback()
  def listen(self, callback):
    if self._listen_handler is None:
      c = self.widget.pressed.connect(self._listener)
      self._listen_handler = c
    self._listeners.append(callback)
  def unlisten(self, callback):
    self._listeners.remove(callback)
  def block(self):
    self._blockedcount += 1
  def unblock(self):
    self._blockedcount -= 1
    if self._blockedcount < 0: self._blockedcount = 0

def find_buttons(widgets,form,toplevel):
  buttons = []
  if "_othertokens" in form._props:
    toks = form._props["_othertokens"]
    for token, args in toks:
      if token == "button":
        widid = str(len(buttons)+1)
        if toplevel: widid = "None-" + widid
        bwidname = ("_button", widid)
        assert bwidname in widgets, (bwidname, widgets.keys())
        bw = buttonwrapper(widgets[bwidname])
        buttons.append(bw)
  return buttons
  
class qtview_primary(object):
  primary = True
  def __init__(self, widget, value_callback, set_callback, listen_callback):
    self.widget = widget
    self._listen_callback = listen_callback
    self._value_callback = value_callback
    self.set = set_callback
    self._listeners = []
    self._listen_handler = None
    self._blockedcount = 0
  def _listener(self, *args):
    #print("LISTEN!", self, self._blockedcount, self._value_callback())
    if self._blockedcount: return
    for callback in self._listeners:
      callback(*args)
  def listen(self, callback):
    if self._listen_handler is None:
      self._listen_handler = self._listen_callback(self._listener)
    self._listeners.append(callback)
  def unlisten(self, callback):
    self._listeners.remove(callback)
  def block(self):
    self._blockedcount += 1
  def unblock(self):
    self._blockedcount -= 1
    if self._blockedcount < 0: self._blockedcount = 0
  def __getattr__(self, attr):
    if attr != "value": raise AttributeError(attr)
    ret = self._value_callback()
    return ret
  def __str__(self):
    return str(self._value_callback())
    
class qtwrap(object):
  def __init__(self, parent, qtui, widgets=None):
    from .anyQt.QtCore import QObject
    assert isinstance(qtui, QObject)
    self.parent = weakref.proxy(parent)
    self.qtui = qtui
    toplevel = False
    if widgets is None: 
      widgets = build_widgets(qtui)
      toplevel = True
    self.widgets = widgets
    self.properties = {}
    self.selectors = {}
    self.is_resources = {}
    form = parent._form    
    formnames = []
    getfunc = None
    if hasattr(form, "type") and form.type in ("none", "hidden"): return
    if form is not None:
      if form.get_membernames() is not None:
        formnames = form.get_membernames()
        getfunc = partial(getattr, form)
      elif form.arraycount > 0:
        if not hasattr(form, "length") or form.length is None: raise ValueError(form.typename)
        formnames = [str(v) for v in range(form.length)]
        getfunc = lambda v: form.__getitem__(int(v))
    self.parent.buttons = find_buttons(widgets, form,toplevel)
     
    for category, name in self.widgets.keys():
      if category == "_widget" and name.find("-") == -1:        
        self.properties[name] = self.widgets[category, name]
        if name not in formnames: raise ValueError(name)
    for category, name in self.widgets.keys():
      if category == "_selector" and name.find("-") == -1:        
        self.selectors[name] = self.widgets[category, name]
        if name not in formnames: raise ValueError(name)
    self.subformfunc = getfunc
    for formname in formnames:
      mform = getfunc(formname)
      if formname in self.properties: 
        if form.arraycount == 0 and mform.is_resource: self.is_resources[formname] = True
        else: self.is_resources[formname] = False
        continue
      view = qtview(mform)
      mwidgets = trim_widgets(widgets, formname+"-")
      view._wrap(qtui, mwidgets)
      if view._wrapped is None:
        if not mform.is_default and not mform._disabled: 
          raise ValueError(formname)
      else:
        self.properties[formname] = view
    parent._do_wrap(self)
    
class qtview(object):
  primary = False
  def __init__(self, obj,form=None):
    import spyder, Spyder
    from Spyder import Object
    if isinstance(obj, Object):
      form = obj._form()
    elif isinstance(obj, spyder.core.spyderform):
      form = obj
    else:
      ok = False
      try: 
        if issubclass(obj, Object):
          form = obj._form()       
          ok = True
      except TypeError:
        pass
      if not ok: raise TypeError(obj)
    
    self._getters = {}
    self._setters = {}
    self._form = form
    self.buttons = []
    self.typename = form.typename
    self._is_resource = False
    if self._form.is_resource: 
      self._is_resource = True 
      self._resourcetype = getattr(Spyder, "Resource" + self.typename + "Array"*form.arraycount)
    self.type = None
    if self.typename is not None and validvar2(self.typename):
      self.type = getattr(Spyder, self.typename + "Array"*form.arraycount)
    self._wrap = partial(qtwrap, self)
    self._listen_callbacks = None
    self._wrapped = None
  def _value(self):
    if self._form.arraycount > 0:
      ret = []
      for n in range(self._form.length):
        v = self._getters.get(str(n), None)
        if v is not None: v = v.value
        ret.append(v)
    else:
      ret = {}
      for prop in self._wrapped.properties:
        ret[prop] = self._getters[prop].value
    return ret    
  def __str__(self):
    return str(self._value())
  def set(self, v):
    if v is not None and self._is_resource:
     vv = self._resourcetype(v).data() 
    elif self.type is None or isinstance(v, self.type): 
      vv = v
    else: 
      vv = self.type(v)
    for prop in self._wrapped.properties:
      if self._form.arraycount > 0:
        try:
          val = vv[int(prop)]
        except IndexError:
          continue
      else:
        val = getattr(vv,prop)
      self._setters[prop](val)
  def _listener(self, dmmy):
    v = self._value()
    for l in self._listen_callbacks:
      l(v)

  def block(self):
    for prop in self._wrapped.properties:
      self._getters[prop].block()

  def unblock(self):
    for prop in self._wrapped.properties:
      self._getters[prop].unblock()
      
  def listen(self, callback):
    if self._listen_callbacks is None:
      self._listen_callbacks = []
      for prop in self._wrapped.properties:
        proplis = self._prop_listeners[prop]
        proplis(self._listener)
    self._listen_callbacks.append(callback)      
  
  def unlisten(self, callback):
    self._listen_callbacks.remove(callback)
    
  def _widget_listener(self, callback, valuefunc, dmmy=None):
    value = valuefunc()
    callback(value)
    
  def _widget_register_listener(self, widget, event, valuefunc, callback):
    f = partial(self._widget_listener, callback, valuefunc)
    signal = getattr(widget, event)
    return signal.connect(f)

  def _optionwidget_listener(self, callback, valuefunc, dmmy):
    value = valuefunc()
    callback(value)

  def _optionwidget_register_listener(self, widget, event, valuefunc, callback):
    f = partial(self._optionwidget_listener, callback, valuefunc)
    signal = getattr(widget, event)
    return signal.connect(f)
  
  def _optionwidget_get(self, options, widget):
    act = widget.currentIndex()
    if act == -1: return None
    return options[act]

  def _optionwidget_set(self, options, widget, value):
    if value is None: 
      widget.setCurrentIndex(-1)
    else:
      assert value in options, (value, options)
      act = options.index(value)
      widget.setCurrentIndex(act)

  def _radiowidget_get(self, options, widgets):
    act = -1
    for wnr, w in enumerate(widgets):
      if w.isChecked(): 
        act = wnr
        break
    if act == -1: return None
    return options[act]

  def _radiowidget_listener(self, w, callback, value):
    if not w.isChecked(): return
    callback(value)

  def _radiowidget_set(self, options, widgets, value):
    if value is None: 
      for w in widgets: 
        widgets[w].setChecked(False)
    assert value in options, (value, options)
    act = options.index(value)
    widgets[act].setChecked(True)

  def _radiowidget_register_listener(self, widgets, event, callback):
    ret = []
    for widget in widgets:
      f = partial(self._radiowidget_listener, widget, callback)
      signal = getattr(widget, event)
      h = signal.connect(f)
      ret.append(h)
    return ret
  
  def _textbuffer_get(self, textbuffer):
    (iter_first, iter_last) = textbuffer.get_bounds()
    return str(textbuffer.get_text(iter_first, iter_last, True))
    
  def _do_wrap(self,wrapped):
    from .anyQt import QtCore
    self._wrapped = wrapped
    self._getters = {}
    self._setters = {}
    self._prop_listeners = {}
    self._wrapwidgets = collect_widgets(wrapped.widgets.values())
    self.widget = find_common_ancestor(self._wrapwidgets)
    def _widget_listener(callback, widget):
      return callback(widget.value)      

    for prop in self._wrapped.properties:
      widget = self._wrapped.properties[prop]      
      wwidget = widget
      mform = self._wrapped.subformfunc(prop)
      resourcewrap = "text"
      if prop in self._wrapped.selectors: resourcewrap = "file"
      if isinstance(widget, qtview):
        self._getters[prop] = widget
        self._setters[prop] = widget.set
        self._prop_listeners[prop] = widget.listen
        if widget.widget is self.widget:
          widget.widget = make_pseudowidget(widget._wrapwidgets, self.widget)
        continue
            
      m_is_resource = self._wrapped.is_resources[prop]
      widgetname = type(widget).__name__
            
      if widgetname in ("QSpinBox", "QDoubleSpinBox"):
        assert not m_is_resource 
        if widgetname == "QSpinBox": typ = int
        elif widgetname == "QDoubleSpinBox": typ = float
        else: raise Exception
        
        getprop = widget.value
        setprop = validval(widget.setValue, typ)
        proplis = partial(self._widget_register_listener, 
         widget, "valueChanged", getprop)
      elif widgetname == "QLineEdit":
        getprop = string_getter(widget.text)
        #setprop = validval(widget.setText, str)
        setprop = partial(lineeditval, widget)        
        if m_is_resource:
          getprop, setprop = wrap_resource(getprop, setprop, mform, resourcewrap)
        proplis = partial(self._widget_register_listener, 
         widget, "textChanged", getprop)
      elif widgetname == "QCheckBox":
        assert not m_is_resource 
        getprop = widget.isChecked
        setprop = validval(widget.setChecked, bool)
        proplis = partial(self._widget_register_listener, 
         widget, "toggled", getprop)
      elif widgetname == "QComboBox":
        assert not m_is_resource 
        options = mform.options
        assert options is not None
        getprop = partial(self._optionwidget_get, options, widget)
        setprop = partial(self._optionwidget_set, options, widget)
        proplis = partial(self._optionwidget_register_listener, 
         widget, "activated", getprop)
      elif widgetname == "QRadioButton":
        assert not m_is_resource 
        options = mform.options
        assert options is not None
        wwidget = [widget]
        for n in range(1,len(options)):
          widg = self._wrapped.widgets["_widget", prop+"-"+str(n)]
          wwidget.append(widg)
        getprop = partial(self._radiowidget_get, options, wwidget)
        setprop = partial(self._radiowidget_set, options, wwidget)
        proplis = partial(self._radiowidget_register_listener, 
         wwidget, "toggled")
      elif widgetname == "QPlainTextEdit":
        getprop = string_getter(widget.toPlainText)
        #setprop = validval(widget.setPlainText, str)
        setprop = partial(texteditval, widget)
        if m_is_resource:
          getprop, setprop = wrap_resource(getprop, setprop, mform, resourcewrap)
        proplis = partial(self._widget_register_listener, 
         widget, "textChanged", getprop)
      elif widgetname == "QTextEdit":
        getprop = string_getter(widget.toPlainText)
        setprop = partial(texteditval, widget)
        if m_is_resource:
          getprop, setprop = wrap_resource(getprop, setprop, mform, resourcewrap)
        proplis = partial(self._widget_register_listener, 
         widget, "textChanged", getprop)
      else:
        n = widget
        try:
          n = widget.get_name()
        except AttributeError:
          pass  
        raise Exception("Widget '%s': Not yet implemented in wrapping" % n)

      wwidgets = trim_widgets(self._wrapped.widgets, prop+"-")
      view = qtview_primary(wwidget,getprop,setprop,proplis)
      view.buttons = find_buttons(wwidgets, mform, False)
      self._getters[prop] = view 
      self._setters[prop] = setprop      
      self._prop_listeners[prop] = proplis

      if prop in self._wrapped.selectors:
        from .anyQt.QtGui import QFileDialog
        selwidget = self._wrapped.selectors[prop]
        if selwidget.property("__file_select__"):
          signal = selwidget.pressed
          dialog_func = partial(openfilename, QFileDialog.getOpenFileName)
          if selwidget.property("__file_select_folder__"):
            dialog_func = QFileDialog.getExistingDirectory
          dialog_args = []
          mformfile = None
          if m_is_resource: 
            mformfile = mform.typename
          elif hasattr(mform, "file"): 
            mformfile = mform.file          
          if mformfile is not None:        
            mformfile += mform.arraycount * "Array"  
            getprop2 = partial(getpropfunc, getprop, mformfile)
            if m_is_resource: getprop2 = getprop #getpropfunc already invoked
            proplis = partial(self._widget_register_listener, 
             widget, "textChanged", getprop2)            
            view._listen_callback = proplis
          slot = partial(selector_eval,dialog_func,dialog_args,setprop)
          signal.connect(slot)
        else:
          raise TypeError("Unknown selector '%s'" % prop)
  
  def __getitem__(self, key):
    assert isinstance(key, int)
    assert key >= 0
    assert self._form.arraycount > 0    
    return self._getters[str(key)]
     
  def __getattr__(self, attr):
    if attr == "value": return self._value()
    try:
      return self._getters[attr]
    except KeyError:
      raise AttributeError(attr)
    
  def __setattr__(self, attr, value):
    if attr.startswith("_") or attr in ("type", "typename","buttons", "widget"): 
      return object.__setattr__(self, attr, value)
    try:
      return self._setters[attr](value)
    except KeyError:
      raise AttributeError(attr)
