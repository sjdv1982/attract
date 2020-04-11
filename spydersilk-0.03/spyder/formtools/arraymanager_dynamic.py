from functools import partial

def add_buttons(form):
  if hasattr(form, "arraymanager") and form.arraymanager == "dynamic":
    assert form.arraycount
    assert hasattr(form, "length")
    for n in range(form.length):
      subform = form[n]
      #KLUDGE: Qtform button adding is still bugged...
      subform = getattr(subform, subform._membernames[0])
      #/KLUDGE
      if n < form.length - 1:
        subform.add_button("Insert element", "before")
        subform.add_button("Append element", "before")
      subform.add_button("Delete element", "before")        
  if form.arraycount:
    if getattr(form,"type",None) in ("none", "hidden"): return
    assert hasattr(form, "length")
    for n in range(form.length):
      subform = form[n]
      add_buttons(subform)
  elif form._membernames is not None:    
    for key in form._membernames:
      subform = getattr(form, key)
      add_buttons(subform)

      
def _insert(controller, index):
  controller._insert(index)
  update_visibility(controller)

def _delete(controller, index):
  controller._delete(index)
  update_visibility(controller)
  
def configure(controller):  
  from .controller import controller_leaf    
  if isinstance(controller, controller_leaf): return
  form = controller._form
  #TODO: hiding stuff!
  if hasattr(form, "arraymanager") and form.arraymanager == "dynamic":
    assert form.arraycount
    assert hasattr(form, "length")
    controller._dynamic_array = True
    for n in range(form.length):
      subcontroller = controller[n]
      subform = form[n]
      #KLUDGE: Qtform button adding is still bugged...
      subcontroller = getattr(subcontroller, subform._membernames[0])
      #/KLUDGE         
      view = subcontroller._view()
      pos = 0
      if n < form.length - 1:
        func1 = partial(_insert, controller, n)
        func2 = partial(_insert, controller, n+1)
        view.buttons[0].listen(func1)
        view.buttons[1].listen(func2)
        pos += 2
      func = partial(_delete, controller, n)
      view.buttons[pos].listen(func)
  if form.arraycount:
    assert hasattr(form, "length")
    for n in range(form.length):
      configure(controller[n])      
  elif form._membernames is not None:    
    for key in form._membernames:
      subcontroller = getattr(controller, key, None)
      if subcontroller is None: continue
      configure(subcontroller)
        
def update_visibility(controller):
  maxchild = 0
  for n in controller._children.keys():
    if n < maxchild: continue
    child = controller._children[n]
    if child._active: maxchild = n 
  form = controller._form
  maxvisible = min(maxchild+1, form.length)
  for n in controller._children.keys():
    child = controller._children[n]
    view = child._view()
    #KLUDGE: Qtform button adding is still bugged...
    child = getattr(child, form[n]._membernames[0])    
    #/KLUDGE             
    view2 = child._view()
    if n < form.length - 1:
      view2.buttons[0].widget.show()
      view2.buttons[1].widget.show()
    if n <= maxvisible:
      view.widget.show()
      if maxvisible == form.length and n < form.length - 1:
        view2.buttons[0].widget.hide()
        view2.buttons[1].widget.hide()  
      elif n == maxvisible and n < form.length - 1:
        view2.buttons[0].widget.hide()
        view2.buttons[1].widget.hide()
        view2.buttons[2].widget.hide()  
      elif n == maxvisible and n == form.length - 1:
        view2.buttons[0].widget.hide()                          
      elif n == maxvisible-1 and n < form.length - 1:
        view2.buttons[1].widget.hide()          
    else:
      view.widget.hide()
    
    
      

def update_all_visibility(controller):
  from .controller import controller_leaf    
  if isinstance(controller, controller_leaf): return  
  form = controller._form
  if hasattr(form, "arraymanager") and form.arraymanager == "dynamic":
    assert form.arraycount
    assert hasattr(form, "length")    
    update_visibility(controller)
    return
  for childname in sorted(controller._children):
    subcontroller = controller._children[childname]
    update_all_visibility(subcontroller)
  