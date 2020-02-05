# Copyright 2008, 2009 Sjoerd de Vries
# This file is part of the Spyder module: "atom" 
# For licensing information, see LICENSE.txt 

import os

naparser = spyder._relative_import("naparser")
whatif = spyder._relative_import("whatif")

ValidationError = spyder.ValidationError
class AtomValidationError(ValidationError):
  pass

