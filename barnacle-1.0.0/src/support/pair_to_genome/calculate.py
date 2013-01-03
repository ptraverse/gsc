#! /usr/bin/env python
"""
calculate.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# TODO
# replace samtools calls with Rod's code?

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
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  CheckConfigCommands)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, EnsureAbsPath, FileBoxCls)
from support.samtools import SAMToolsCls
from group import P2GGroupCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "P2G SUPPORT SUCCESS"
MSG_FAIL = "P2G SUPPORT FAIL"

class P2GCalculatorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    CheckConfigCommands(self, "samtools")
    self.groups_file = None
    self.output_file = None
    self.options.use_chr = False
  #} end def

  def __del__(self): #{
    # close input and output files, if they are not already closed
    self.CloseFiles()
    CloseLogFile(self)
  #} end def

  def CalculateSupport(self): #{
    start = time.time()
    LogMsg(self, "Adding pair-to-genome support to groups...")
    # open the input and output files
    self.Setup()
    #ExtremeDebugMsg(self, "Should I use chr? %s" % self.options.use_chr)
    # for each group in the input file
    for group_line in self.groups_file: #{
      group_start = time.time()
      # create a group object from the line
      group = P2GGroupCls(self.options, self.log_info)
      group.ParseGroupLine(group_line)
      LogMsg(self, "Group: %i" % group.group_id)
      ExtremeDebugMsg(self, "  %s" % group.ToString())
      # get the pair-to-genome support for the current group
      group.GetPairToGenomeSupport()
      # write the pair-to-genome support for the current group
      self.WritePairToGenomeSupport(group.SupportString())
      ExtremeDebugMsg(self, "Time spent on group: %s" % TimeSpent(group_start))
    #} end for
    # close the input and output files
    self.CloseFiles()
    # remove the temporary samtools output files
    for end in ["", "_1", "_2"]: #{
      temp_sam_path = os.path.join(self.options.output_dir,
        "sam_out_tmp%s" % end)
      if (os.path.isfile(temp_sam_path)): #{
        os.remove(temp_sam_path)
      #} end if
      temp_sam_path += ".err"
      if (os.path.isfile(temp_sam_path)): #{
        os.remove(temp_sam_path)
      #} end if
    #} end for
    LogMsg(self, "Total time adding pair-to-genome support: %s" %
      TimeSpent(start))
  #} end def

  def Setup(self): #{
    fail_msg = "cannot open groups file"
    self.groups_file = FileBoxCls(self.options.barnacle_path, "r", fail_msg)
    output_file_path = self.options.barnacle_path.replace(".data", ".out")
    fail_msg = "cannot create pair-to-genome support output file"
    self.output_file = FileBoxCls(output_file_path, "w", fail_msg)
    # create samtools object and check whether to use "chr" in chromosome IDs
    samtools = SAMToolsCls(self.options.p2g_path, self.options,
      log_info=self.log_info)
    self.options.use_chr = samtools.ShouldChromUseChr()
  #} end def

  def WritePairToGenomeSupport(self, support_string): #{
    self.output_file.WriteLine("%s" % support_string)
  #} end def

  def CloseFiles(self): #{
    if (None != self.groups_file and not self.groups_file.closed): #{
      self.groups_file.close()
      self.groups_file = None
    #} end if
    if (None != self.output_file and not self.output_file.closed): #{
      self.output_file.close()
      self.output_file = None
    #} end if
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class P2GCalculatorError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Calculates pair-to-genome support for groups in "
    "the input file")
  args = [ "BARNACLE_FILE", "PAIR_TO_GENOME_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--read-length",
                    type="int", metavar="N",
                    help="Each read is N bp long. [ default: %default ]")
  parser.add_option("--frag-len",
                    type="int",
                    help="The expected length of the mate-pair fragments "
                         "is N bp. [ default: %default ]")
  parser.add_option("--frag-fract",
                    type="float", metavar="F",
                    help="When looking for read-pairs spanning the genomic "
                         "event region, only count a pair as being "
                         "significantly different from the expected length if "
                         "the fractional difference between the observed and "
                         "expected fragment lengths is more than F. "
                         "[ default: %default ]")
  parser.add_option("--min-mapq",
                    type="int", metavar="N",
                    help="When looking for read-pairs spanning the genomic "
                         "event region, filter pairs with mapping quality "
                         "less than N. [ default: %default ]")
  parser.add_option("--max-sam-retries",
                    type="int", metavar="N",
                    help="If there is an error running a samtools view "
                         "command, retry the command a maximum of N times. "
                         "[ default: %default ]")
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("--log-file",
                    dest="log_file_name", metavar="FILE",
                    help="Log all messages in FILE")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.set_defaults(read_length = 50,
                      frag_len    = 200,
                      frag_fract  = 0.10,
                      min_mapq    = 10,
                      max_sam_retries = 5,
                      dpt=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  path_errors = list()
  CheckFilePath(options.barnacle_path, "groups", path_errors)
  # use the directory of the groups file for the output directory
  options.output_dir = os.path.dirname(options.barnacle_path)
  CheckFilePath(options.p2g_path, "pair-to-genome alignments", path_errors)
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n%s" % "\n".join(path_errors))
  #} end if
  # the paths are good if there are no path errors
  return (0 == len(path_errors))
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.barnacle_path = EnsureAbsPath(args[0])
    options.p2g_path      = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try:
        p2g_calculator = P2GCalculatorCls(options)
        WriteCommand(p2g_calculator, sys.argv)
        p2g_calculator.CalculateSupport()
      except (MyError), e:
        ErrMsg("ERROR while calculating pair-to-genome support:\n  %s" % e)
        return ES_RUN_ERR
      # end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify the path to a Barnacle pair-to-genome "
                 "data file (BARNACLE_FILE); and the path to a bam-formatted "
                 "read-pair to genome alignment file (PAIR_TO_GENOME_FILE).")
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
