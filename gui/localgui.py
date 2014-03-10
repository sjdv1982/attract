from functools import partial
import attractsave


def buttons(form, obj):
  form.add_button("Generate shell script", "after")
  form.add_button("Save", "after")

  
def logic(m,v,c, silk):
  v.buttons[0].listen(partial(silk.save, m))
  v.buttons[1].listen(partial(attractsave.generate, m, silk.outputfile))
  v.listen(silk.viewupdate)
  silk.viewupdate()  
  
  
  