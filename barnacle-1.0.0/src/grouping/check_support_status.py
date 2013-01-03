#! /usr/bin/env python
"""
check_support_status.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, signal, sys, time, traceback

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
from utils.log import CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  GetOptions)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  FileBoxCls, CleanLine)
from utils.subprocesses import (RunCommandFromString, RunPipeFromStrings,
  STRING_OUT)
from utils.library_iterator import LibIteratorCls
from alignment_processing import check_status as check_cid_status
from support.pair_to_genome import check_status as check_p2g_status
from support.read_to_contig import check_status as check_r2c_status

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "CHECK STATUS SUCCESS"
MSG_FAIL = "CHECK STATUS FAIL"

class CheckSupportStatusCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    self.status_checkers = dict()
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def CheckAllLibraries(self): #{
    DebugMsg(self, "Checking status of libraries...")
    old_handler = signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    lib_list_file = FileBoxCls(self.options.lib_list_path, "r",
      "could not open library list file")
    # do not check event paths while iterating
    self.options.no_path_check = True
    # get the default options for each library
    self.lib_options = {
      'cid': GetOptions(check_cid_status),
      'p2g': GetOptions(check_p2g_status),
      'r2c': GetOptions(check_r2c_status),
    }
    self.lib_options['p2g'].terse = True
    self.lib_options['r2c'].terse = True
    lib_iterator = LibIteratorCls(self.options.lib_list_path,
       self.CheckLibraryStatus, self.options, self.log_info)
    LogMsg(self, "="*80)
    lib_iterator.IterateOverAllLibs()
    LogMsg(self, "Checked support status for %i libraries" %
      lib_iterator.num_libs)
    signal.signal(signal.SIGPIPE, old_handler)
  #} end def

  #  for lib_line in lib_list_file: #{
  #    lib_line = lib_line.rstrip("\n\r\t ")
  #    # skip blank or comment lines
  #    if (0 == len(lib_line) or lib_line.startswith("#")): #{
  #      continue
  #    #} end if
  #    self.status = None
  #    LogMsg(self, "="*80)
  #    self.GetLogFilePath(lib_line)
  #    LogMsg(self, self.lib_name)
  #    if (self.options.check_cid): #{
  #      self.CheckCIDStatus()
  #    #} end if
  #    self.CheckLibSupportStatus()
  #    self.CheckP2GStatus()
  #    self.DisplayStatus()
  #  #} end for
  #  lib_list_file.close()
  #  LogMsg(self, "="*80)
  #  signal.signal(signal.SIGPIPE, old_handler)
  #} end def

  def CheckLibraryStatus(self, lib_info): #{
    self.status = None
    LogMsg(self, "Checking status of %s..." % lib_info.lib_name)
    support_log_path = self.GetLogFilePath(lib_info)
    if (self.options.check_cid): #{
      self.CheckCIDStatus(lib_info.lib_dir)
    #} end if
    if (None != support_log_path): #{
      self.CheckSupportStatus(support_log_path)
    #} end if
    if (self.options.check_p2g): #{
      self.CheckP2GStatus(lib_info.lib_dir)
    #} end if
    if (self.options.check_r2c): #{
      self.CheckR2CStatus(lib_info.lib_dir)
    #} end if
    self.DisplayStatus(support_log_path, lib_info.lib_dir)
    LogMsg(self, "="*80)
  #} end def

  def GetLogFilePath(self, lib_info): #{
    #(self.lib_name, self.assembly_ver, self.barnacle_ver) = lib_line.split(",")
    #self.support_dir = os.path.join(
    #  self.options.libs_dir,
    #  self.lib_name,
    #  "Assembly",
    #  "assembly-%s" % self.assembly_ver,
    #  "anomalous_contigs",
    #  "ver_%s" % self.barnacle_ver,
    #)
    #self.support_dir = self.options.lib_path_template
    #self.support_dir = self.support_dir.replace("DIR_HEAD",
    #  self.options.libs_dir)
    #self.support_dir = self.support_dir.replace("%{lib}", self.lib_name)
    #if ("current" != self.assembly_ver and
    #    not self.assembly_ver.startswith("assembly-")):
    #  self.assembly_ver = "assembly-" + self.assembly_ver
    ##} end if
    #if ("current" != self.barnacle_ver and
    #    not self.barnacle_ver.startswith("ver_")):
    #  self.barnacle_ver = "ver_" + self.barnacle_ver
    ##} end if
    #self.support_dir = self.support_dir.replace("%{assembly_ver}", self.assembly_ver)
    #self.support_dir = self.support_dir.replace("%{barnacle_ver}", self.barnacle_ver)
    #self.support_log_file_name = "%s.anomalous_contig.log.continue" % self.lib_name
    #support_log_path = os.path.join(self.support_dir, "1_raw_candidates",
    #  support_log_file_name)
    support_log_file_name = "%s.barnacle.log.continue" % lib_info.lib_name
    support_log_path = os.path.join(lib_info.lib_dir, "1_raw_candidates",
      support_log_file_name)
    DebugMsg(self, "Support log file path: %s" % support_log_path)
    if (not os.path.isfile(support_log_path)): #{
      self.status = "not started / no log file"
      #raise CheckSupportStatusError \
      #  ("cannot find log file: %s" % support_log_path)
      return None
    #} end if
    return support_log_path
  #} end def

  def CheckSupportStatus(self, support_log_path): #{
    if (None != self.status): #{
      return
    #} end if
    DebugMsg(self, "Checking support status...")
    tail_cmd = "tail -n1 %s" % support_log_path
    output = CleanLine(RunCommandFromString(tail_cmd, stdout=STRING_OUT,
      dpt=self.options.dpt)[0])
    #DebugMsg(self, "STATUS: \"%s\"" % output)
    if ("complete" == output.lower()): #{
      self.status = "complete"
    else:
      self.status = "in progress"
    #} end if
  #} end def

  def CheckCIDStatus(self, lib_dir): #{
    start = time.time()
    #LogMsg(self, "Checking candidate identification status...")
    self.lib_options['cid'].jobs_dir = os.path.join(lib_dir, "cluster_cid")
    try: #{
      cid_checker = check_cid_status.CIDStatusCls(self.lib_options['cid'],
        log_info=self.log_info)
      cid_checker.CheckCIDStatus()
      self.status_checkers['cid'] = cid_checker
    except MyError, e:
      LogMsg(self, "Error checking candidate identification status: %s" % e)
    #} end try
    LogMsg(self, "Time spent checking candidate identification status: %s" %
      TimeSpent(start))
  #} end def

  def CheckP2GStatus(self, lib_dir): #{
    start = time.time()
    #LogMsg(self, "Checking pair-to-genome status...")
    #p2g_options = self.options
    #p2g_options.debug = False
    p2g_jobs_dir = os.path.join(lib_dir, "cluster_p2g")
    self.lib_options['p2g'].jobs_dir = p2g_jobs_dir
    if (not os.path.exists(p2g_jobs_dir)): #{
      LogMsg(self, "PAIR-TO-GENOME...\n  not parallel OR not yet submitted")
      DebugMsg(self, "  %s" % p2g_jobs_dir)
      return
    #} end if
    #p2g_options.resubmit = False
    try: #{
      p2g_checker = check_p2g_status.P2GStatusCheckerCls(
        self.lib_options['p2g'], log_info=self.log_info)
      p2g_checker.CheckP2GStatus()
      self.status_checkers['p2g'] = p2g_checker
    except MyError, e:
      LogMsg(self, "Error checking pair-to-genome status: %s" % e)
    #} end try
    #p2g_checker.OutputP2GStatus()
    LogMsg(self, "Time spent checking pair-to-genome status: %s" %
      TimeSpent(start))
  #} end def

  def CheckR2CStatus(self, lib_dir): #{
    start = time.time()
    #LogMsg(self, "Checking read-to-contig status...")
    r2c_jobs_dir = os.path.join(lib_dir, "cluster_r2c")
    self.lib_options['r2c'].jobs_dir = r2c_jobs_dir
    if (not os.path.exists(r2c_jobs_dir)): #{
      LogMsg(self, "READ-TO-CONTIG...\n  not parallel OR not yet submitted")
      DebugMsg(self, "  %s" % r2c_jobs_dir)
      return
    #} end if
    try: #{
      r2c_checker = check_r2c_status.R2CStatusCheckerCls(
        self.lib_options['r2c'], log_info=self.log_info)
      r2c_checker.CheckR2CStatus()
      self.status_checkers['r2c'] = r2c_checker
    except MyError, e:
      LogMsg(self, "Error checking read-to-contig status: %s" % e)
    #} end try
    LogMsg(self, "Time spent checking read-to-contig status: %s" %
      TimeSpent(start))
  #} end def

  def DisplayStatus(self, support_log_path, lib_dir): #{
    if (sys.stdout.closed): #{
      DebugMsg(self, "stdout closed!")
      return
    #} end if
    DebugMsg(self, "Barnacle status...")
    #old_handler = signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    #LogMsg(self, "="*80)
    if ('cid' in self.status_checkers and
        None != self.status_checkers['cid']): #{
      LogMsg(self, "CANDIDIATE IDENTIFICATION...")
      self.status_checkers['cid'].OutputCIDStatus()
      self.status_checkers['cid'] = None
    #} end if
    if ('p2g' in self.status_checkers and
        None != self.status_checkers['p2g']): #{
      LogMsg(self, "PAIR-TO-GENOME...")
      self.status_checkers['p2g'].OutputP2GStatus()
      self.status_checkers['p2g'] = None
    #} end if
    if ('r2c' in self.status_checkers and
        None != self.status_checkers['r2c']): #{
      LogMsg(self, "READ-TO-CONTIG...")
      self.status_checkers['r2c'].OutputR2CStatus()
      self.status_checkers['r2c'] = None
    #} end if
    LogMsg(self, "SUPPORT STATUS: %s" % self.status)
    if ("complete" == self.status): #{
      LogMsg(self, "LOCATION: %s" % lib_dir)
    elif ("in progress" == self.status):
      pattern = r'"\.\.\.$"'
      grep_cmd = "grep %s %s" % (pattern, support_log_path)
      tail_cmd = "tail -n1"
      status = CleanLine(RunPipeFromStrings([grep_cmd, tail_cmd],
        stdout=STRING_OUT, dpt=self.options.dpt)[0])
      LogMsg(self, status)
      #grep_process = subprocess.Popen(
      #    shlex.split(grep_cmd),
      #    stdout=subprocess.PIPE,
      #    preexec_fn=lambda:
      #      signal.signal(signal.SIGPIPE, signal.SIG_DFL))
      #status = subprocess.Popen(
      #    shlex.split(tail_cmd),
      #    stdin=grep_process.stdout,
      #    stdout=subprocess.PIPE).communicate()[0]
      #LogMsg(self, status, newline=False)
    #} end if
    sys.stdout.flush()
    #signal.signal(signal.SIGPIPE, old_handler)
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class CheckSupportStatusError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Checks the status of Barnacle support jobs submitted "
    "to the cluster.")
  args = [ "LIB_LIST_FILE", "DIR_HEAD", ]
  usage_string = ("%prog " + " ".join(args) + " [ OPTIONS ]\n"
    "LIB_LIST_FILE is a text file containing a list of library info lines, "
    "each line of the form \"lib_name,assembly_ver,barnacle_ver\" (e.g. "
    "\"A00001,assembler-1.2.1,1.1.0\").\nDIR_HEAD is the directory containing "
    "the subdirectories of the libraries in LIB_LIST_FILE.")
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-c", "--check-cid",
                    action="store_true", dest="check_cid",
                    help="Check the status of the candidate identification "
                         "jobs as well as the support job.")
  parser.add_option("-p", "--check-p2g",
                    action="store_true", dest="check_p2g",
                    help="Check the status of the pair-to-genome support "
                         "jobs as well as the support job.")
  parser.add_option("-r", "--check-r2c",
                    action="store_true", dest="check_r2c",
                    help="Check the status of the read-to-contig support "
                         "jobs as well as the support job.")
  parser.add_option("-t", "--template",
                    metavar="TEMPLATE", dest="lib_path_template",
                    help="The template to use to construct the path to the "
                         "event data files for a library. \"DIR_HEAD\", "
                         "\"%{lib}\", \"%{assembly_ver}\", and "
                         "\"%{barnacle_ver}\" will be replaced with the "
                         "appropriate values from LIB_LIST_FILE. "
                         "[default: %default]")
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("--log-file",
                    dest="log_file_name", metavar="FILE",
                    help="Log all messages in FILE")
  parser.add_option("--terse",
                    action="store_true", dest="terse",
                    help="When checking candidate identification, "
                         "pair-to-genome, or read-to-contig status, only "
                         "report number complete and total number of jobs.")
  parser.add_option("-q", "--quiet",
                    action="store_true",
                    help="Only write output to log-file, not to the screen.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(check_cid=False,
                      check_p2g=False,
                      check_r2c=False,
                      lib_path_template= os.path.join(["DIR_HEAD", "%{lib}",
                        "Assembly", "%{assembly_ver}", "barnacle",
                        "%{barnacle_ver}"]),
                      dpt=False,
                      terse=False,
                      quiet=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  path_errors = list()
  CheckFilePath(options.lib_list_path, "library list", path_errors)
  CheckDirPath(options.libs_dir, "library base", path_errors)
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n%s" % "\n".join(path_errors))
  #} end if
  # require log file in quiet mode
  if (options.quiet and None == options.log_file_name): #{
    ErrMsg("please specify a log file to use when using the quiet option")
    opts_good = False
  #} end if
  # the paths are good if there are no path errors
  return (opts_good and 0 == len(path_errors))
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.lib_list_path = EnsureAbsPath(args[0])
    options.libs_dir      = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try: #{
        status_checker = CheckSupportStatusCls(options)
        WriteCommand(status_checker, sys.argv)
        status_checker.CheckAllLibraries()
      except MyError, e:
        ErrMsg("ERROR while checking status:\n  %s" % e)
        return ES_RUN_ERR
      except IOError, e:
        ErrMsg("Caught IOError:\n  %s" % e)
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
