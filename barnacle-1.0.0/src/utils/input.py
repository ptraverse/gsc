#! /usr/bin/env python
"""
input.py

Created by Lucas Swanson 2012 09 19
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import readline, sys

# import custom modules
from utils.error import MyError
from utils.messages import LogMsg

def GetStringInput(self, prompt="", write_to_log=True): #{
  #print "Getting string input!"
  sys.stdout.flush()
  sys.stderr.flush()
  #if (hasattr(self, "log_info")): #{
  #  self.log_info["file"].Flush()
  #} end if
  raw_value = raw_input(prompt)
  #print "Got value: %s" % raw_value
  if (write_to_log): #{
    LogMsg(self, "%s%s" % (prompt, raw_value), write_to_screen=False)
  #} end if
  return raw_value
#} end def

def GetYesOrNoInput(self, prompt=""): #{
  #print "Getting yes or no input!"
  raw_value = GetStringInput(self, prompt, write_to_log=False)
  value = raw_value[0].upper()
  #print "Value converted to: %s" % value
  while (value not in "YN"): #{
    raw_value = GetStringInput(self, "  Please enter \"yes\" or \"no\": ",
        write_to_log=False)
    value = raw_value[0].upper()
  #} end while
  LogMsg(self, "%s%s" % (prompt, value), write_to_screen=False)
  return value
#} end def

#### EXCEPTION CLASSES ####
class UserInputError(MyError): #{
  """Exception raised for errors encountered while getting user input"""
  pass
#} end class
