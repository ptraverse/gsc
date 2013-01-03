#! /usr/bin/env python
"""
multi_lib_check_status.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, sys, traceback

# import custom modules
from version import VERSION
from utils.log import CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, WriteCommand, GetOptions)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath)
from utils.library_iterator import LibIteratorCls
import check_status

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

class MultiLibCIDStatusCls:
  def __init__(self, options): #{
    SetupMainClass(self, options)
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def CheckAllLibraries(self): #{
    # do not check event paths while iterating
    self.options.no_path_check = True
    # get the default options for each library
    self.lib_options = GetOptions(check_status)
    lib_iterator = LibIteratorCls(self.options.lib_list_path,
       self.CheckLibraryStatus, self.options, self.log_info)
    lib_iterator.IterateOverAllLibs()
    LogMsg(self, "Checked CID status for %i libraries" %
      lib_iterator.num_libs)
    #fail_msg = "could not open library list file"
    #lib_list_file = OpenFile(lib_list_path, "r", fail_msg)
    #for lib_line in lib_list_file: #{
    #  lib_line = lib_line.rstrip("\n\r ")
    #  if (0 == len(lib_line) or lib_line.startswith("#")): #{
    #    continue
    #  #} end if
    #  (lib, dir) = lib_line.split(",")
    #  print lib
    #  if (not os.path.isdir(dir)): #{
    #    raise MultiLibCIDStatusError("\"%s\" is not a directory" % dir)
    #  #} end if
    #  status_checker = check_status.CIDStatusCls(self.options)
    #  try: #{
    #    status_checker.CheckStatus(dir)
    #    status_checker.Output()
    #  except MyError, e:
    #    log_msg(self.log_file, "Error checking library status: %s" % e)
    #  #} end try
    #} end for
    #lib_list_file.close()
  #} end def

  def CheckLibraryStatus(self, lib_info): #{
    self.lib_options.jobs_dir = os.path.join(lib_info.lib_dir, "cluster_cid")
    try: #{
      LogMsg(self, lib_info.lib_name)
      status_checker = check_status.CIDStatusCls(self.lib_options,
        log_info=self.log_info)
      status_checker.CheckCIDStatus()
      status_checker.OutputCIDStatus()
      LogMsg(self, "")
    except MyError, e:
      LogMsg(self, "Error checking library status: %s" % e)
    #} end try
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class MultiLibCIDStatusError(MyError):
  pass
#} end class

def SetupOptionsParser():
  description_string = ("Runs the check ACF status script on all libraries "
    "in the input file.")
  args = [ "LIB_LIST_FILE", "DIR_HEAD", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--log-file",
                    dest="log_file_name", metavar="FILE",
                    help="Log all messages in FILE")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  path_errors = list()
  CheckFilePath(options.lib_list_path, "library list", path_errors)
  CheckDirPath(options.libs_dir, "library base", path_errors)
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n%s" % "\n".join(path_errors))
  #} end if
  # the paths are good if there are no path errors
  return (0 == len(path_errors))# and opts_good)
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if(parser.num_args == len(args)): #{
    options.lib_list_path = EnsureAbsPath(args[0])
    options.libs_dir      = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try: #{
        full_status_checker = MultiLibCIDStatusCls(options)
        WriteCommand(full_status_checker, sys.argv)
        full_status_checker.CheckAllLibraries()
      except MyError, e:
        ErrMsg("ERROR while checking status:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library list file (LIB_LIST_FILE); "
      "and a base directory to use when constructing the library paths "
      "(DIR_HEAD).")
    return ES_OPT_ERR
  #} end if
  return ES_SUCCESS
#} end def

if __name__ == '__main__':
  try: #{
    exit_status = Main()
  except Exception, e:
    traceback.print_exc()
    exit_status = ES_EXCEPTION
  #} end try
  if (ES_SUCCESS == exit_status): #{
    print MSG_SUCCESS
  else:
    print MSG_FAIL
  #} end if
  sys.exit(exit_status)
#} end if
