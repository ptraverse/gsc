#! /usr/bin/env python
"""
submit.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, re, sys, time, traceback

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
from general import (SetupMainClass, TimeSpent, WriteCommand,
  CheckConfigCommands, ConfigHasCommand, GetCommand, GetClusterValue)
from messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from files_paths import CheckFilePath, EnsureAbsPath
from subprocesses import RunCommandFromString

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "JOB SUBMISSION SUCCESS"
MSG_FAIL = "JOB SUBMISSION FAIL"

class JobSubmitterCls: #{
  def __init__(self, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    if (not options.single): #{
      self.mqsub_cmd = os.path.join(os.environ["BARNACLE_PATH"],
        "utils", "mqsub")
      CheckFilePath(self.mqsub_cmd, "mqsub binary")
    #  CheckConfigCommands(self, "mqsub")
    #  #CheckConfigCommands(self, "mqsub", local=False, check_path=False,
    #  #  cluster_head=self.options.cluster_head)
    #} end if
    if (None == self.options.hostname): #{
      self.options.hostname = GetClusterValue(self, "hostname")
    #} end if
    if (None == self.options.queue): #{
      self.options.queue = GetClusterValue(self, "queue")
    #} end if
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def Submit(self): #{
    LogMsg(self, "Submitting job(s)...")
    start = time.time()
    qsub_args = self.GetQSubArgs()
    if (self.options.single): #{
      command = self.SingleJobCommand(qsub_args)
    else:
      command = self.MultipleJobsCommand(qsub_args)
    #} end if
    self.SubmitCommand(command)
    LogMsg(self, "Time spent submitting job(s): %s" % TimeSpent(start))
  #} end def

  def GetQSubArgs(self): #{
    qsub_args_list = [
      "-l hostname=%s" % self.options.hostname,
      "-q", self.options.queue,
      "-o", self.options.job_dir,
      "-e", self.options.job_dir,
    ]
    for mem_field in ["mem_free", "mem_token", "h_vmem"]: #{
      mem_arg = "-l %s=%s" % (mem_field, self.options.mem)
      qsub_args_list.append(mem_arg)
    #} end for
    if (None != self.options.wall_time): #{
      qsub_args_list.append("-l walltime=%s" % self.options.wall_time)
    #} end if
    if (None != self.options.email): #{
      qsub_args_list.append("-m a -M %s" % self.options.email)
    #} end if
    return " ".join(qsub_args_list)
  #} end def

  def SingleJobCommand(self, qsub_args): #{
    command_list = [
      "qsub",
      "-N", self.options.job_name,
      qsub_args,
      "-wd", self.options.job_dir,
      self.options.job_path,
    ]
    return " ".join(command_list)
  #} end def

  def MultipleJobsCommand(self, qsub_args): #{
    #mqsub_cmd = GetCommand(self, "mqsub")
    #mqsub_cmd = GetCommand(self, "mqsub", local=False,
    #  cluster_head=self.options.cluster_head)
    command_list = [
      self.mqsub_cmd,
      "--name", self.options.job_name,
      "--qsub \'%s\'" % qsub_args,
      "--chdir", self.options.job_dir,
      "--file", self.options.job_path,
    ]
    return " ".join(command_list)
  #} end def

  def SubmitCommand(self, command): #{
    #if (ConfigHasCommand(self, "cluster_setup")): #{
    #  command = ("%s && %s" % (GetCommand(self, "cluster_setup"), command))
    #} end if
    cluster_setup = GetClusterValue(self, "setup")
    if (cluster_setup.lower() not in ["", "none"]): #{
      command = "source %s && %s" % (cluster_setup, command)
    #} end if
    submit_command = "ssh %s \"%s\"" % (self.options.cluster_head, command)
    DebugMsg(self, "Submit Command:\n  %s" % submit_command)
    returncode = RunCommandFromString(submit_command, dpt=self.options.dpt)
    if (0 != returncode): #{
      raise JobSubmitterError("error submitting job, return code: %i\n%s" %
        (returncode, submit_command))
    #} end if
  #} end def
#} end class

def SetSubmitOptions(submit_options, job_name, path,
    single, from_options): #{
  submit_options.job_name     = job_name
  submit_options.job_path     = path
  submit_options.single       = single
  submit_options.cluster_head = from_options.cluster_head
  submit_options.mem          = from_options.mem
  submit_options.hostname     = from_options.hostname
  submit_options.queue        = from_options.queue
  submit_options.wall_time    = from_options.wall_time
  submit_options.debug        = from_options.debug
#} end def

def GetJobNum(job_dir): #{
  job_match = re.search(r"[^0-9](?P<num>\d+)$", job_dir)
  if (None == job_match): #{
    raise JobSubmitterError("cannot get job number from "
      "job directory %s" % job_dir)
  #} end if
  return int(job_match.group('num'))
#} end def

#### EXCEPTION CLASSES ####
class JobSubmitterError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = "Submit jobs to cluster"
  args = [ "JOB_NAME", "JOB_FILE", "CLUSTER_HEAD", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-s", "--single",
                    action="store_true",
                    help="JOB_FILE is a single job [default].")
  parser.add_option("-m", "--multiple",
                    action="store_false", dest="single",
                    help="JOB_FILE is a list of jobs.")
  parser.add_option("--mem", "--memory",
                    help="The memory requirement of the jobs "
                         "[default:\"%default\"].")
  parser.add_option("--hostname",
                    help="The hostname(s) to submit to "
                         "[default: read from barnacle.cfg].")
                         #"[default:\"%default\"].")
  parser.add_option("--queue",
                    help="The queue(s) to submit to [default: read from "
                      "barnacle.cfg].")
                    #help="The queue(s) to submit to [default:\"%default\"].")
  parser.add_option("-w", "--wall-time",
                    metavar="H:MM:SS",
                    help="The maximum time to spend on the job.")
  parser.add_option("--email",
                    help="E-mail status updates to the given email address")
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(single=True,
                      mem="1G",
                      dpt=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  # check the memory
  mem_pattern = r"^(?P<value>\d+)(?P<unit>[GgMmKk])$"
  mem_match = re.search(mem_pattern, options.mem)
  if (None == mem_match): #{
    ErrMsg("  Invalid memory option: \"%s\". Must be integer value followed "
      "by one of G, M, or K. e.g. \"1G\"")
    opts_good = False
  else:
    try:
      options.mem = "%i%s" % \
        (int(mem_match.group('value')), mem_match.group('unit'))
    except ValueError, e:
      ErrMsg("  Invalid memory option: \"%s\". Must be integer value followed "
        "by one of G, M, or K. e.g. \"1G\"")
      opts_good = False
    # end try
  #} end if
  bad_chars = [r"\n", r"\t", r"\r", r"/", r":", r"@", "\\", r"*", r"?"]
  if (None != re.search(r"^[0-9]", options.job_name) or
      None != re.search(r"[\n\t\r/:@\*?]", options.job_name)): #{
    ErrMsg("  Invalid job name \"%s\": cannot start with " % options.job_name +
      "a digit or contain any of: %s." % ", ".join(["\"%s\"" % bad_char for
      bad_char in bad_chars]))
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.job_path, "job(s)", path_errors)
  options.job_dir = os.path.dirname(options.job_path)
  options.log_file_name = GetLogPath(options.job_path,
    "submit.%s" % options.job_name, options.job_dir)
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # the paths are good if there are no path errors and no conflicting options
  return (opts_good and 0 == len(path_errors))
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.job_name     = args[0]
    options.job_path     = EnsureAbsPath(args[1])
    options.cluster_head = args[2]
    if (CheckPaths(options)): #{
      try:
        submitter = JobSubmitterCls(options)
        WriteCommand(submitter, sys.argv)
        submitter.Submit()
      except (MyError), e:
        ErrMsg("ERROR while submitting job(s):\n  %s" % e)
        return ES_RUN_ERR
      # end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a job name (JOB_NAME), the path to a "
      "job(s) file (JOB_FILE), and a cluster (CLUSTER_HEAD) to submit to.")
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
