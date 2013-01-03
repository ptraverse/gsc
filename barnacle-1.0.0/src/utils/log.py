#! /usr/bin/env python
"""
log.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os

# import custom modules
from messages import DebugMsg
from files_paths import FileBoxCls
from subprocesses import STRING_OUT, RunCommandFromString

def GetLogPath(input_path, task_desc, output_dir=None): #{
  (input_dir, input_file_name) = os.path.split(input_path)
  if (None == output_dir): #{
    output_dir = input_dir
  #} end if
  (input_root, ext) = os.path.splitext(input_file_name)
  # set the log-file name
  if ("" == task_desc): #{
    log_file_name = "%s.log" % (input_root)
  else:
    log_file_name = "%s.%s.log" % (input_root, task_desc)
  #} end if
  return os.path.join(output_dir, log_file_name)
#} end def

def OpenLogFile(log_file_name, mode="a"): #{
  if (None == log_file_name): #{
    return None
  #} end if
  fail_msg = "Could not create log file"
  log_file = FileBoxCls(log_file_name, mode, fail_msg)
  log_file.Write("START LOG: %s\n" % RunCommandFromString \
    ("date", stdout=STRING_OUT, dpt=True)[0])
  return log_file
#} end def

def CloseLogFile(self): #{
  if (self.close_log_file): #{
    #DebugMsg(self, "Closing log file...")
    if (None != self.log_file): #{
      self.log_file.close()
      self.log_file = None
    #} end if
  #} end if
#} end def
