# Copyright 2011, Sjoerd de Vries
# This file is part of the Spyder module: "core" 
# For licensing information, see LICENSE.txt 

class spyderlist(list,Object):
  """Wrapper class around a Python list
  The conversion engine does not operate on List"""    
  @staticmethod
  def typename(): return "List"  
  def __init__(self, *a):
    try:
      list.__init__(self,*a)
    except:
      list.__init__(self,a)
    self.__validate__()
  def __validate__(self):
    pass      
  def dict(self): 
    """Called by the dict function of Spyder classes
    that have List members
    For internal use only"""
    return self
  def __print__(self, spaces, mode):
    """Called by the pretty-print function of Spyder classes
    that have List members
    For internal use only"""        
    return str(self)
