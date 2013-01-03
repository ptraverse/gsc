#! /usr/bin/env python
"""
tokenizer.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
#from utils.messages import DebugMsg

# CONSTANTS

class TokenizerCls: #{
  def __init__(self, input_string, delimiter="\t", log_info=None): #{
    self.log_info = log_info
    self.delimiter = delimiter
    self.fields = input_string.strip(delimiter).split(delimiter)
    self.index = 0
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  def next(self): #{
    #DebugMsg(self, "Getting next token")
    if (len(self.fields) <= self.index): #{
      raise StopIteration
    #} end if
    self.index += 1
    #DebugMsg(self, "  Token: %s" % self.fields[self.index-1])
    return self.fields[self.index-1]
  #} end def

  def Reset(self): #{
    self.index = 0
  #} end def

  def Skip(self, num=1):
    self.index += num
    #DebugMsg(self, "Skipping %i tokens" % num)
    #for i in range(num): #{
    #  token = self.next()
    #  DebugMsg(self, "Skipping token: %s" % token)
    #} end def
  #} end def

  def Size(self): #{
    return len(self.fields)
  #} end def

  def InputString(self): #{
    return self.delimiter.join(self.fields)
  #} end def
#} end class

class ParseWarningCls: #{
  def __init__(self, warning, data, rule=""): #{
    self.warning = warning
    self.data    = data
    self.rule    = rule
  #} end def
#} end class

def GetFieldAndValue(field_str): #{
  return field_str.split(":", 1)
#} end def

def GetFieldValue(field_str): #{
  try:
    return GetFieldAndValue(field_str)[1]
  except IndexError, e:
    raise TokenizerError \
      ("Error getting field value from \"%s\": %s" % (field_str, e))
  # end try
#} end def

#### EXCEPTION CLASSES ####
class TokenizerError(MyError): #{
  pass
#} end class
