#! /usr/bin/env python
"""
check_format.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, sys, time, traceback

# ensure that the sys.path is set up to allow importing Barnacle modules
barnacle_dir = sys.path[0]
if ("" == barnacle_dir): #{
  barnacle_dir = os.getcwd()
#} end if
while (not os.path.isfile(os.path.join(barnacle_dir, "barnacle.pl"))): #{
  barnacle_dir = os.path.dirname(barnacle_dir)
  if ("" == barnacle_dir or "/" == barnacle_dir): #{
    print "cannot find BARNACLE directory"
    raise Exception
  #} end if
#} end while
barnacle_dir = os.path.expanduser(barnacle_dir)
barnacle_dir = os.path.abspath(barnacle_dir)
if (barnacle_dir not in sys.path): #{
  sys.path.insert(0, barnacle_dir)
#} end if

# import custom modules
from version import VERSION
from log import GetLogPath, CloseLogFile
from error import MyError
from general import SetupMainClass, TimeSpent, WriteCommand
from messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, GetFilePath, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "CHECK EVENTS SUCCESS"
MSG_FAIL = "CHECK EVENTS FAIL"

class EventCheckerCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    self.output_files = dict()
    self.num_bad_events = 0
  #} end def

  def __del__(self): #{
    for output_file in self.output_files.itervalues(): #{
      output_file.close()
    #} end for
    CloseLogFile(self)
  #} end def

  def CheckForUnparsableEvents(self): #{
    LogMsg(self, "Removing events with errors...")
    start = time.time()
    # create a file to hold the unparsable events
    self.CreateOutputFiles()
    # create an event parser
    parser = CandidateGroupParserCls(self.options.barnacle_path,
      keep_lines=True, check_data=True)
    # loop through the events
    for event in parser: #{
      DebugMsg(self, "Event %i, Warnings: %i, Member Warnings: %i" %
        (event.id, len(event.warnings), len(event.member_warnings)))
      # if there were no parsing warnings
      if (0 == (len(event.warnings) + len(event.member_warnings))): #{
        # write the event to the good output
        self.output_files['good'].Write(event.FullDataString())
      else:
        DebugMsg(self, event.WarningsString())
        DebugMsg(self, "\n".join(event.member_warnings))
        # write the event lines to the unparsable output
        self.output_files['bad'].WriteLine("%s" % parser.GroupLine())
        if (0 < len(event.warnings)): #{
          self.output_files['bad'].WriteLine("PARSING ERRORS:\n%s" %
            event.WarningsString(indent="  "))
        #} end if
        for (member_line, member) in zip (parser.MemberLines(),
            event.members): #{
          self.output_files['bad'].WriteLine("%s" % member_line)
          if (0 < len(member.warnings)): #{
            self.output_files['bad'].WriteLine("PARSING ERRORS:\n%s" %
              member.WarningsString(indent="  "))
          #} end if
        #} end for
        # write an extra line between events
        self.output_files['bad'].WriteLine("")
        # increment the number of unparsable events
        self.num_bad_events += 1
      #} end if
    #} end for
    # output warnings about unparsable events
    self.WarnAboutUnparsableEvents()
    LogMsg(self, "Time spent removing events with errors: %s" %
      TimeSpent(start))
  #} end def

  def CreateOutputFiles(self): #{
    good_events_path = GetFilePath(self.options.output_dir,
      self.options.lib, "data")
    DebugMsg(self, "Parsable events path: \"%s\"" % good_events_path)
    fail_msg = "cannot create parsable events file"
    self.output_files['good'] = FileBoxCls(good_events_path, "w", fail_msg)
    bad_events_path = GetFilePath(self.options.output_dir,
      self.options.lib, "unparsable")
    DebugMsg(self, "Unparsable events path: \"%s\"" % bad_events_path)
    fail_msg = "cannot create unparsable events file"
    self.output_files['bad'] = FileBoxCls(bad_events_path, "w", fail_msg)
  #} end def

  def WarnAboutUnparsableEvents(self): #{
    if (0 == self.num_bad_events): #{
      LogMsg(self, "No unparsable events found!")
      return
    #} end if
    # if errors, message: file containing events with errors
    LogMsg(self, "WARNING: %i unparsable events " % self.num_bad_events +
      "found. They have been moved to: %s" % self.output_files['bad'].path)
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class EventCheckerError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Makes a backup of the input file, and creates a new "
    "file with any unparsable events removed. Unparsable events are saved in "
    "their own file.")
  args = [ "LIB", "BARNACLE_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  # get and check the output path
  options.output_dir = GetOutDir(os.path.dirname(options.barnacle_path),
    "parsable_candidates")
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors, create=True)
    # get the log-file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "check_format", options.output_dir)
  #} end if
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # the paths are good if there are no path errors and no conflicting options
  return (0 == len(path_errors))
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.lib           = args[0]
    options.barnacle_path = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try:
        event_checker = EventCheckerCls(options)
        WriteCommand(event_checker, sys.argv)
        event_checker.CheckForUnparsableEvents()
      except (MyError), e:
        ErrMsg("ERROR while checking for unparsable events:\n  %s" % e)
        return ES_RUN_ERR
      # end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); and the path to a "
      "Barnacle data file for that library (BARNACLE_FILE).")
    return ES_OPT_ERR
  #} end if
  return ES_SUCCESS
#} end def

if __name__ == '__main__': #{
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
