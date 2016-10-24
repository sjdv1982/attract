import sys

class TVLog(object):
  print_log = False
  def __init__(self):
    self.reset()
  def reset(self):
    self.plan = None
    self.text = ""
  def set_plan(self, plan):
    self.plan = plan
    if self.print_log: 
      print >> sys.stderr, self.plan
  def log(self, *args):
    line = " ".join([str(a) for a in args])
    self.text += line + "\n"
    if self.print_log: 
      print >> sys.stderr, line
    