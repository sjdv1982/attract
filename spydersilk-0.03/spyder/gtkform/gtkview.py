from functools import partial
import weakref

def if_validval(func, typ, value):
  if value == None: return
  val = typ(value)
  return func(val)

def validval(func, typ):
  return partial(if_validval, func, typ)

def getpropfunc(getprop, filetype):
  return (getprop(), filetype)

def filegetstr(widget):
  fn = widget.true_filename
  if fn is None: return None
  assert isinstance(fn, str), fn
  return fn
  
def fileset(widget, value):  
  widget.unselect_all()
  if value is None or value == "": return    
  widget.true_filename = value
  widget.set_filename(value)  
  
def file_update(widget, *args, **kwargs):
  widget.true_filename = widget.get_filename()
  
def wrap_resource(getprop, setprop, form, formtype):
  from spyder.formtools.resourcewrapper import ResourceWrapper
  typename = form.typename + form.arraycount * "Array"
  typename2 = "Resource"+typename
  getprop2 = getprop
  if formtype == "file": 
    getprop2 = partial(getpropfunc, getprop, typename)
  w = ResourceWrapper(getprop2, setprop, typename2, formtype)
  return w.getprop, w.setprop
  
def build_widgets(gtkbuilder):
  from gi.repository import Gtk  
  ret = {}
  widgets = gtkbuilder.get_objects()
  for w in widgets:
    try:
      name = Gtk.Buildable.get_name(w)
    except TypeError:
      continue
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
      c = self.widget.connect("button-press-event", self._listener)
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
  
class gtkview_primary(object):
  primary = True
  def __init__(self, widget, value_callback, set_callback, listen_callback):
    self.widget = widget
    self._listen_callback = listen_callback
    self._value_callback = value_callback
    self._set_callback = set_callback
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
  def set(self, value):
    #Gtk widgets annoyingly delay the update of their contents:
    # widget.set_value("x), widget.get_value() typically does not work
    # Therefore, we have to block our listeners when we update manually
    self.block()
    self._set_callback(value)
    self.unblock()
  def block(self):
    self._blockedcount += 1
  def unblock(self):
    self._blockedcount -= 1
    if self._blockedcount < 0: self._blockedcount = 0
  def __getattr__(self, attr):
    if attr != "value": raise AttributeError(attr)
    return self._value_callback()
  def __str__(self):
    return str(self._value_callback())
    
class gtkwrap(object):
  def __init__(self, parent, gtkbuilder, widgets=None):
    from gi.repository import Gtk
    assert isinstance(gtkbuilder, Gtk.Builder)
    self.parent = weakref.proxy(parent)
    self.gtkbuilder = gtkbuilder
    toplevel = False
    if widgets is None: 
      widgets = build_widgets(gtkbuilder)
      toplevel = True
    self.widgets = widgets
    self.properties = {}
    self.is_resources = {}
    form = parent._form    
    formnames = []
    getfunc = None
    if hasattr(form, "type") and form.type == "none": return
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
    self.subformfunc = getfunc
    for formname in formnames:
      mform = getfunc(formname)
      if formname in self.properties: 
        if form.arraycount == 0 and mform.is_resource: self.is_resources[formname] = True
        else: self.is_resources[formname] = False
        continue
      view = gtkview(mform)
      mwidgets = trim_widgets(widgets, formname+"-")
      view._wrap(gtkbuilder, mwidgets)
      if view._wrapped is None:
        if not mform.is_default: raise ValueError(formname)
      else:
        self.properties[formname] = view
    if len(self.properties):
      parent._do_wrap(self)
    
class gtkview(object):
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
    self.type = None
    if self.typename is not None:
      self.type = getattr(Spyder, self.typename + "Array"*form.arraycount)
    self._wrap = partial(gtkwrap, self)
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
    if self.type is None: raise AttributeError
    if v is not None and self._is_resource:
     vv = self._resourcetype(v).data() 
    elif isinstance(v, self.type): 
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
    
  def _widget_listener(self, callback, valuefunc, w):
    value = valuefunc()
    callback(value)

  def _file_widget_register_listener(self, widget, filetype):
    def getprop(*args, **kwargs):
      return widget.true_filename, filetype
    proplis = partial(self._widget_register_listener, 
     widget, "file-set", getprop)
    return proplis 
    
  def _widget_register_listener(self, widget, event, valuefunc, callback):
    f = partial(self._widget_listener, callback, valuefunc)
    return widget.connect(event, f)

  def _optionwidget_listener(self, callback, valuefunc, dmmy):
    value = valuefunc()
    callback(value)

  def _optionwidget_register_listener(self, widget, event, valuefunc, callback):
    f = partial(self._optionwidget_listener, callback, valuefunc)
    return widget.connect(event, f)
  
  def _optionwidget_get(self, options, widget):
    act = widget.get_active()
    if act == -1: return None
    return options[act]

  def _optionwidget_set(self, options, widget, value):
    if value is None: 
      widget.set_active(-1)
    else:
      assert value in options, (value, options)
      act = options.index(value)
      widget.set_active(act)

  def _radiowidget_get(self, options, widgets):
    act = -1
    for wnr, w in enumerate(widgets):
      if w.get_active(): 
        act = wnr
        break
    if act == -1: return None
    return options[act]

  def _radiowidget_listener(self, callback, valuefunc, w):
    if not w.get_active(): return
    value = valuefunc()
    callback(value)

  def _radiowidget_set(self, options, widgets, value):
    if value is None: 
      for w in widgets: 
         widgets[w].set_active(False)
    else:    
      assert value in options, (value, options)
      act = options.index(value)
      widgets[act].set_active(True)

  def _radiowidget_register_listener(self, widgets, event, valuefunc, callback):
    ret = []
    for widget in widgets:
      f = partial(self._radiowidget_listener, callback, valuefunc)
      h = widget.connect(event, f)
      ret.append(h)
    return ret
  
  def _textbuffer_get(self, textbuffer):
    (iter_first, iter_last) = textbuffer.get_bounds()
    return str(textbuffer.get_text(iter_first, iter_last, True))
    
  def _do_wrap(self,wrapped):
    self._wrapped = wrapped
    self._getters = {}
    self._setters = {}
    self._prop_listeners = {}
    from gi.repository import Gtk  
    def _widget_listener(callback, widget):
      return callback(widget.value)      

    for prop in self._wrapped.properties:
      widget = self._wrapped.properties[prop]
      wwidget = widget
      mform = self._wrapped.subformfunc(prop)      
      if isinstance(widget, gtkview):
        self._getters[prop] = widget
        self._setters[prop] = widget.set
        self._prop_listeners[prop] = widget.listen
        continue
      
      m_is_resource = self._wrapped.is_resources[prop]
      
      try:
        widgetname = widget.get_name()
      except AttributeError:
        widgetname = None
        
      if widgetname == "GtkSpinButton":
        assert not m_is_resource 
        getprop = widget.get_value
        setprop = validval(widget.set_value, float)
        proplis = partial(self._widget_register_listener, 
         widget, "value-changed", getprop)
      elif widgetname == "GtkEntry":
        if m_is_resource:
          getprop, setprop = wrap_resource(getprop, setprop, mform, "text")        
        getprop = widget.get_text
        setprop = validval(widget.set_text, str)
        proplis = partial(self._widget_register_listener, 
         widget, "changed", getprop)
      elif widgetname == "GtkFileChooserButton":
        widget.true_filename = None
        getprop = partial(getattr, widget, "true_filename")
        widget.connect("file-set", partial(file_update, widget)) 
        setprop = partial(fileset, widget)
        
        mformfile = None
        if m_is_resource: 
          mformfile = mform.typename
        elif hasattr(mform, "file"): 
          mformfile = mform.file          
        if mformfile is not None:        
          mformfile += mform.arraycount * "Array"  
        
        if mformfile is not None:        
          if m_is_resource:
            getprop2 = partial(filegetstr, widget)
            getprop, setprop = wrap_resource(getprop2, setprop, mform, "file")
          else:  
            getprop = partial(getpropfunc, getprop, mformfile)
          proplis = self._file_widget_register_listener(widget, mformfile)          
        else:
          proplis = partial(self._widget_register_listener, 
           widget, "file-set", getprop)
      elif widgetname == "GtkCheckButton":
        assert not m_is_resource 
        getprop = widget.get_active
        setprop = validval(widget.set_active, int)
        proplis = partial(self._widget_register_listener, 
         widget, "toggled", getprop)
      elif widgetname == "GtkComboBoxText":
        assert not m_is_resource 
        options = mform.options
        assert options is not None
        getprop = partial(self._optionwidget_get, options, widget)
        setprop = partial(self._optionwidget_set, options, widget)
        proplis = partial(self._optionwidget_register_listener, 
         widget, "changed", getprop)
      elif widgetname == "GtkRadioButton":
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
         wwidget, "toggled", getprop)
      elif widgetname == "GtkTextView":
        buf = widget.get_buffer()
        getprop = partial(self._textbuffer_get, buf)
        setprop = validval(buf.set_text, str)
        if m_is_resource:
          getprop, setprop = wrap_resource(getprop, setprop, mform, "text")        
        proplis = partial(self._widget_register_listener, 
         buf, "changed", getprop)
      else:
        n = widget
        try:
          n = widget.get_name()
        except AttributeError:
          pass  
        raise Exception("Widget '%s': Not yet implemented in wrapping" % n)

      wwidgets = trim_widgets(self._wrapped.widgets, prop+"-")
      view = gtkview_primary(wwidget,getprop,setprop,proplis)
      view.buttons = find_buttons(wwidgets, mform, False)
      self._getters[prop] = view 
      self._setters[prop] = setprop      
      self._prop_listeners[prop] = proplis
  
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
    if attr.startswith("_") or attr in ("type", "typename","buttons"): 
      return object.__setattr__(self, attr, value)
    try:
      return self._setters[attr](value)
    except KeyError:
      raise AttributeError(attr)
