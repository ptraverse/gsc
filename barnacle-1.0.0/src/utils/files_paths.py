#! /usr/bin/env python
"""
files_paths.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import errno, os, re

#import custom modules
from error import MyError

def GetFilePath(dir, lib, ext):
  file_name = "%s.barnacle.%s" % (lib, ext)
  return os.path.join(dir, file_name)
# end def

def CheckFilePath(file_path, file_desc, path_errors=None): #{
  if (not os.path.isfile(file_path)): #{
    msg = "Cannot find %s file: \"%s\"." % (file_desc, file_path)
    if (None == path_errors): #{
      raise PathError(msg)
    else:
      path_errors.append(msg)
    #} end if
  #} end if
#} end def

def CheckNewFilePath(file_path, file_desc, path_errors=None): #{
  if (os.path.exists(file_path)): #{
    msg = "%s file already exists: \"%s\"." % (file_desc, file_path)
    if (None == path_errors): #{
      raise PathError(msg)
    else:
      path_errors.append(msg)
    #} end if
  elif (not os.path.isdir(os.path.dirname(file_path))): #{
    msg = "%s directory does not exist: \"%s\"." % (file_desc,
      os.path.dirname(file_path))
    if (None == path_errors): #{
      raise PathError(msg)
    else:
      path_errors.append(msg)
    #} end if
  #} end if
#} end def

def CheckDirPath(dir_path, dir_desc, path_errors=None,
    create=False, replace=True): #{
  if (create and not replace and os.path.exists(dir_path)): #{
    msg = ("The %s directory already exists, please move the old directory "
      "so that it will not be overwritten or use the --force option to "
      "force overwriting:\n%s." % (dir_desc, dir_path))
    if (None == path_errors): #{
      raise PathError(msg)
    else:
      path_errors.append(msg)
      return
    #} end if
  #} end if
  task = "find"
  if (create): #{
    task = "create"
    EnsureDirectoryExists(dir_path)
  #} end if
  if (not os.path.isdir(dir_path)): #{
    msg = "Cannot %s %s directory: \"%s\"." % (task, dir_desc, dir_path)
    if (None == path_errors): #{
      raise PathError(msg)
    else:
      path_errors.append(msg)
    #} end if
  #} end if
#} end def

def EnsureDirectoryExists(path): #{
  if (os.path.exists(path)): #{
    return
  #} end if
  try: #{
    os.makedirs(path)
  except OSError,e:
    if (errno.EEXIST != e.errno): #{
      raise e
    #} end if
  #} end try
#} end def

def EnsureAbsPath(path): #{
  path = os.path.expanduser(path)
  if (os.path.isabs(path)): #{
    return path
  #} end if
  return os.path.abspath(path)
#} end def

def GetOutDir(in_dir, out_tail): #{
  clean_in_dir = in_dir.rstrip("/")
  in_base = os.path.basename(clean_in_dir)
  in_match = re.search(r"^(\d+)(\.\d+|[a-z])*_.*", in_base)
  if (None == in_match): #{
    out_num = 1
    out_parent = clean_in_dir
  else:
    out_num = int(in_match.group(1)) + 1
    out_parent = os.path.dirname(clean_in_dir)
  #} end if
  out_dir = os.path.join(out_parent, "%i_%s" % (out_num, out_tail))
  return out_dir
#} end def

def OpenFile(filename, mode, fail_msg="could not open file"): #{
  try:
    file = open(filename, mode)
  except IOError, e:
    raise FileError(fail_msg, filename, e)
  # end try
  return file
#} end def

def CleanLine(line): #{
  return line.rstrip("\n\r\t ")
#} end def

class FileBoxCls: #{
  def __init__(self,
      file_path, open_mode, open_fail_message="could not open file",
      clean_lines=True, skip_blank_lines=True):
    self.path = file_path
    self.mode = open_mode
    self.open_fail_message = open_fail_message
    self.clean_lines = clean_lines
    self.skip_blank_lines = skip_blank_lines
    self.file = None
    self.closed = True
  #} end def

  def __del__(self): #{
    self.Close()
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  def Open(self): #{
    self.Close()
    self.file = OpenFile(self.path, self.mode, self.open_fail_message)
    self.closed = self.file.closed
  #} end def

  def Write(self, to_write): #{
    if (None == self.file or self.file.closed): #{
      self.Open()
    #} end if
    self.file.write(to_write)
  #} end def

  def WriteLine(self, to_write): #{
    self.Write("%s\n" % to_write)
  #} end def

  def next(self): #{
    if (None == self.file or self.file.closed): #{
      self.Open()
    #} end if
    line = self.file.next()
    if (self.clean_lines): #{
      line = CleanLine(line)
    #} end if
    if (self.skip_blank_lines and "" == line): #{
      return self.next()
    #} end if
    return line
  #} end def

  def Close(self): #{
    if (None != self.file and not self.file.closed): #{
      self.file.close()
      self.closed = self.file.closed
    #} end if
  #} end def

  def close(self): #{
    self.Close()
  #} end def

  def Head(self, num_lines=1): #{
    head_lines = list()
    for line in self: #{
      head_lines.append(line)
      if (len(head_lines) >= num_lines): #{
        break
      #} end if
    #} end for
    self.Close()
    if (1 == num_lines): #{
      return head_lines[0]
    #} end if
    return head_lines[:num_lines]
  #} end def

  def Tail(self, num_lines=1): #{
    tail_lines = list()
    for i in range(num_lines): #{
      tail_lines.append(None)
    #} end for
    for line in self: #{
      for i in range(len(tail_lines)-1): #{
        tail_lines[i] = tail_lines[i+1]
      #} end for
      tail_lines[-1] = line
    #} end for
    self.Close()
    if (1 == num_lines): #{
      return line
    #} end if
    return tail_lines
  #} end def

  def Flush(self): #{
    if (None != self.file and not self.closed): #{
      self.file.flush()
    #} end if
  #} end def
#} end class

class FileError(MyError): #{
  def __init__(self, fail_msg, filename, error): #{
    msg = "%s: %s\n    %s" % (fail_msg, filename, error)
    MyError.__init__(self, msg)
  #} end def
#} end class

class PathError(MyError): #{
  pass
#} end class
