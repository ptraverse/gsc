#! /usr/bin/env python
"""
messages.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import sys

# import custom modules
from files_paths import FileBoxCls

# constants
TEST_MODE = False

def ErrMsg(msg, newline=True): #{
  """Write an error message to stderr (appends a newline)"""
  if (newline): #{
    print >>sys.stderr, msg
  else:
    print >>sys.stderr, msg,
  #} end if
#} end def

def TestMsg(msg, newline=True): #{
  if (TEST_MODE): #{
    ErrMsg(msg, newline)
  #} end if
#} end if

def LogMsg(target, msg, newline=True, write_to_screen=True): #{
  if (None != target and
      hasattr(target, 'options') and
      hasattr(target.options, 'quiet') and
      target.options.quiet): #{
    write_to_screen = False
  #} end if
  if (write_to_screen): #{
    ErrMsg(msg, newline)
  #} end if
  if (None == target): #{
    return
  #} end if
  log_file = None
  if (isinstance(target, FileBoxCls)): #{
    log_file = target
  if (isinstance(target, dict) and 'file' in target):
    log_file = target['file']
  elif (hasattr(target, 'log_info') and
      None != target.log_info and
      'file' in target.log_info):
    log_file = target.log_info['file']
  elif (hasattr(target, 'log_file')):
    log_file = target.log_file
  #} end if
  if (None != log_file): #{
    log_file.Write(msg)
    if (newline): #{
      log_file.Write("\n")
    #} end if
  #} end if
#} end def

def DebugMsg(target, msg, newline=True): #{
  if (None == target): #{
    return
  #} end if
  debug = False
  if (isinstance(target, dict) and 'debug' in target): #{
    debug = target['debug']
  elif (hasattr(target, "log_info") and
      None != target.log_info and
      'debug' in target.log_info):
    debug = target.log_info['debug']
  elif (hasattr(target, "debug")):
    debug = target.debug
  #} end if
  if (debug): #{
    LogMsg(target, msg, newline)
  #} end if
#} end def

def ExtremeDebugMsg(target, msg, newline=True): #{
  if (None == target): #{
    return
  #} end if
  debug = False
  if (isinstance(target, dict) and 'extreme_debug' in target): #{
    debug = target['extreme_debug']
  elif (hasattr(target, "log_info") and
      None != target.log_info and
      'extreme_debug' in target.log_info):
    debug = target.log_info['extreme_debug']
  elif (hasattr(target, "extreme_debug")):
    debug = target.extreme_debug
  #} end if
  if (debug): #{
    LogMsg(target, msg, newline)
  #} end if
#} end def
