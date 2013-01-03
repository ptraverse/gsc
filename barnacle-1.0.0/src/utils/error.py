#! /usr/bin/env python
"""
error.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

class MyError(Exception): #{
  def __init__(self, msg): #{
    self.msg = msg
  #} end def
  def __str__(self): #{
    return str(self.msg)
  #} end def
#} end class

class ParsingError(MyError): #{
  pass
#} end class
