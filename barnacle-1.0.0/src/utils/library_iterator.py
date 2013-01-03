#! /usr/bin/env python
"""
library_iterator.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os

# import custom modules
from utils.error import MyError
from utils.messages import DebugMsg
from utils.files_paths import FileBoxCls

class LibIteratorCls:
  def __init__(self, lib_list_path, ProcessLibraryMethod,
      options, log_info=None):
    self.lib_list_path = lib_list_path
    self.ProcessLibrary = ProcessLibraryMethod
    self.options = options
    self.CheckOptions()
    #if (not hasattr(options, "no_path_check")):
    #  self.options.no_path_check = False
    # end if
    self.log_info = log_info
    self.num_libs = 0
    self.list_of_paths = False
  # end def

  def __del__(self):
    if (hasattr(self, "lib_list_file") and
        None != self.lib_list_file     and
        not self.lib_list_file.closed):
      self.lib_list_file.close()
    # end if
  # end def

  def CheckOptions(self): #{
    required_opts = ["no_path_check", "only_pass"]
    for opt in required_opts: #{
      if (not hasattr(self.options, opt)): #{
        setattr(self.options, opt, False)
      #} end if
    #} end for
    if (not hasattr(self.options, "get_paths")):
      self.options.get_paths = True
    #} end if
  #} end def

  def IterateOverAllLibs(self):
    self.num_libs = 0
    self.lib_list_file = FileBoxCls(self.lib_list_path, "r",
      "could not open library list file")
    for lib_line in self.lib_list_file:
      # skip comment lines
      if (lib_line.startswith("#")):
        continue
      # end if
      lib_info = LibraryInfoCls(self.options, self.log_info,
        self.list_of_paths)
      lib_info.GetLibDir(lib_line)
      DebugMsg(self, "Lib Dir: %s" % lib_info.lib_dir)
      if (self.options.get_paths and not self.list_of_paths):
        lib_info.GetEventPaths()
      # end if
      if (1 > len(lib_info.event_paths) and
          (self.options.get_paths or self.list_of_paths)):
        raise LibIteratorError("could not get event path(s) from directory: "
          "%s" % lib_info.lib_dir)
      # end if
      self.ProcessLibrary(lib_info)
      self.num_libs += 1
    # end for
    self.lib_list_file.close()
  # end def
# end class

class LibraryInfoCls:
  def __init__(self, options, log_info=None, list_of_paths=False):
    self.lib_name = None
    self.assembly_ver = None
    self.barnacle_ver = None
    self.lib_dir = None
    self.options = options
    self.log_info = log_info
    self.use_filtered_results = False
    self.event_paths = list()
    self.list_of_paths = list_of_paths
  # end def

  def GetLibDir(self, lib_line):
    DebugMsg(self, "LIB LINE: \"%s\"" % lib_line)
    lib_fields = lib_line.split(",")
    (self.lib_name, self.assembly_ver, self.barnacle_ver) = lib_fields[:3]
    if (self.list_of_paths):
      self.lib_dir = os.path.dirname(lib_fields[3])
      self.event_paths.append(lib_fields[3])
      return
    # end if
    if (None == self.options.libs_dir):
      self.lib_dir = self.options.lib_path_template
    else:
      self.lib_dir = os.path.join(self.options.libs_dir,
        "%{lib}/Assembly/%{assembly_ver}/barnacle/%{barnacle_ver}")
    # end if
    self.lib_dir = self.lib_dir.replace("%{lib}", self.lib_name)
    self.lib_dir = self.lib_dir.replace("%{assembly_ver}", self.assembly_ver)
    if ("current" != self.barnacle_ver and
        not self.barnacle_ver.startswith("ver_")):
      self.barnacle_ver = "ver_" + self.barnacle_ver
    # end if
    self.lib_dir = self.lib_dir.replace("%{barnacle_ver}", self.barnacle_ver)
    # add the filter directory, if any is found
    if (3 < len(lib_fields)):
      self.lib_dir = os.path.join(self.lib_dir, lib_fields[3])
      self.use_filtered_results = True
    # end if
    if (not self.options.no_path_check and
        not os.path.isdir(self.lib_dir)):
      raise LibIteratorError("could not find library directory: %s" %
        self.lib_dir)
    # end if
  # end def

  def GetEventPaths(self):
    if (self.options.only_pass or self.use_filtered_results):
      passing_path = os.path.join(self.lib_dir,
        "%s.barnacle.pass" % self.lib_name)
      self.event_paths.append(passing_path)
      if (not self.options.only_pass):
        fail_path = os.path.join(self.lib_dir,
          "%s.barnacle.fail" % self.lib_name)
        self.event_paths.append(fail_path)
        self.options.append = True
      # end if
    else:
      data_path = os.path.join(self.lib_dir,
        "%s.barnacle.data" % self.lib_name)
      self.event_paths.append(data_path)
    # end if
    for event_path in self.event_paths:
      if (not self.options.no_path_check and
          not os.path.isfile(event_path)):
        raise LibIteratorError("could not find event file: %s" % event_path)
      # end if
    # end for
  # end def
# end class

#### EXCEPTION CLASSES ####
class LibIteratorError(MyError):
  pass
#} end class
