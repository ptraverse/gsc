#! /usr/bin/env python
"""
check_status.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
from glob import iglob
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
from utils import submit
from utils.log import CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  GetOptions, GetClusterValue)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckDirPath, EnsureAbsPath, CheckNewFilePath,
  FileBoxCls)
from utils.subprocesses import RunCommandFromString
from calculate import MSG_SUCCESS as P2G_SUCCESS
from calculate import MSG_FAIL as P2G_FAIL

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "CHECK STATUS SUCCESS"
MSG_FAIL = "CHECK STATUS FAIL"

class P2GStatusCheckerCls: #{
  def __init__(self, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    self.ResetCounts()
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

  def ResetCounts(self): #{
    self.num_jobs = {
      'not started': 0,
      'in progress': 0,
      'failed':      0,
      'complete':    0,
      'total':       0,
    }
    self.unfinished_jobs = list()
  #} end def

  def CheckP2GStatus(self): #{
    LogMsg(self, "Checking pair-to-genome jobs status...\n%s" %
      self.options.jobs_dir)
    self.ResetCounts()
    start_time = time.time()
    job_paths = self.GetJobPaths()
    # for each job in the jobs directory
    for job_path in job_paths: #{
      job = P2GJobCls(job_path)
      if (job.NotReallyJob()): #{
        continue
      #} end if
      self.num_jobs['total'] += 1
      DebugMsg(self, "Job file: %s" % job_path)
      if (not job.GetOutputPath()): #{
        self.num_jobs['not started'] += 1
        DebugMsg(self, "  - not started")
        self.unfinished_jobs.append(job)
        continue
      #} end if
      DebugMsg(self, "Checking job %s, output file: %s" %
        (job.num, job.output_path))
      job.CheckStatus()
      if ("in progress" != job.status): #{
        self.num_jobs[job.status] += 1
      #} end if
      if ("complete" != job.status): #{
        self.unfinished_jobs.append(job)
      #} end if
      DebugMsg(self, "  - %s" % job.status)
    #} end for
    if (0 == self.num_jobs['total']): #{
      raise P2GStatusCheckerError("could not find any jobs in \"%s\"" %
        self.options.jobs_dir)
    #} end if
    self.num_jobs['in progress'] = self.num_jobs['total'] - (
      self.num_jobs['failed'] +
      self.num_jobs['complete'] +
      self.num_jobs['not started'])
    self.SetStatus()
    if (not self.options.terse): #{
      LogMsg(self, "Time spent checking status: %s" % TimeSpent(start_time))
    #} end if
  #} end def

  def GetJobPaths(self): #{
    # construct the job path template
    job_path_template = os.path.join(self.options.jobs_dir, "*-p2g.*.sh")
    job_paths = iglob(job_path_template)
    return job_paths
  #} end def

  def SetStatus(self): #{
    if (0 < self.num_jobs['in progress']): #{
      self.status = "in progress"
    elif (0 < self.num_jobs['failed']): #{
      self.status = "failed"
    elif (0 < self.num_jobs['not started']): #{
      self.status = "waiting"
    elif (self.num_jobs['total'] == self.num_jobs['complete']): #{
      self.status = "complete"
    else:
      self.status = "unknown"
    #} end if
  #} end def

  def ResubmitJobs(self): #{
    LogMsg(self, "Resubmitting unfinished jobs...")
    # check that the right number of jobs were put into
    # the unfinished jobs list
    num_unfinished_jobs = \
      self.num_jobs['in progress'] + \
      self.num_jobs['failed'] + \
      self.num_jobs['not started']
    if (len(self.unfinished_jobs) != num_unfinished_jobs): #{
      raise P2GStatusCheckerError("wrong number of unfinished jobs in list")
    #} end if
    for unfinished_job in self.unfinished_jobs: #{
      self.ResubmitJob(unfinished_job)
    #} end for
    LogMsg(self, "All unfinished jobs resubmitted")
  #} end def

  def ResubmitJob(self, job): #{
    if (None == self.options.cluster_head): #{
      self.ResubmitLocal(job)
    else:
      self.ResubmitCluster(job)
    #} end if
  #} end def

  def ResubmitLocal(self, job): #{
    resubmit_cmd = job.path
    fail_msg = "cannot create stdout file for job resubmission"
    out_file = {'path':job.NewOutPath(), 'mode':"w", 'fail_msg':fail_msg}
    fail_msg = "cannot create stderr file for job resubmission"
    err_file = {'path':job.NewErrPath(), 'mode':"w", 'fail_msg':fail_msg}
    DebugMsg(self, resubmit_cmd)
    return_code = \
      RunCommandFromString(resubmit_cmd, stdout=out_file, stderr=err_file,
        dpt=self.options.dpt)
    if (0 != return_code and not self.log_info['debug']): #{
      LogMsg(self, resubmit_cmd)
    #} end if
    if (0 > return_code): #{
      raise P2GStatusCheckerError("Resubmission command was terminated "
        "by signal %i" % return_code)
    elif (0 < return_code):
      raise P2GStatusCheckerError("Error running resubmission "
        "command: %i" % return_code)
    #} end if
  #} end def

  def ResubmitCluster(self, job): #{
    submit_options = GetOptions(submit)
    submit.SetSubmitOptions(submit_options, job.JobName(), job.path,
      True, self.options)
    if (submit.CheckPaths(submit_options)): #{
      submitter = submit.JobSubmitterCls(submit_options,
        log_info=self.log_info)
      submitter.Submit()
    else:
      raise CheckCandidateIdentificationStatusError("cannot resubmit jobs: "
        "errors in resubmit options")
    #} end if
  #} end def

  def OutputP2GStatus(self): #{
    if (self.options.terse): #{
      status_list = [
        'complete',
        'total'
      ]
    else:
      status_list = [
        'not started',
        'in progress',
        'failed',
        'complete',
        'total'
      ]
    #} end if
    for status_type in status_list: #{
      LogMsg(self, "%s: %i" % (status_type, self.num_jobs[status_type]))
    #} end for
    LogMsg(self, "P2G status: %s" % self.status)
  #} end def
#} end class

class P2GJobCls: #{
  def __init__(self, job_path): #{
    self.path = job_path
    file_name = os.path.basename(job_path)
    # ensure that the path is really a job path
    job_match = re.search(r"^(?P<lib>.*)-p2g\.(?P<num>\d+)\.sh$", file_name)
    if (None == job_match): #{
      self.num = -1
    else:
      self.lib = job_match.group('lib')
      self.num = job_match.group('num')
    #} end if
  #} end def

  def NotReallyJob(self): #{
    if (not hasattr(self, 'num') or None == self.num or -1 == self.num): #{
      return True
    else:
      return False
    #} end if
  #} end def

  def GetOutputPath(self): #{
    self.output_path = None
    output_paths = iglob("%s.o*" % self.path)
    for path in output_paths: #{
      job_id_pattern = r"o(?P<id>\d+)$"
      job_id_match = re.search(job_id_pattern, path)
      if (None == job_id_match): #{
        raise P2GStatusCheckerError(
          "cannot get job id from job output file: %s" % path)
      #} end if
      job_id = int(job_id_match.group("id"))
      # get the most recent output file for the job
      # i.e. the path with the largest job id number
      if (None == self.output_path or job_id > max_job_id): #{
        self.output_path = path
        max_job_id = job_id
      #} end if
    #} end for
    if (None == self.output_path): #{
      return False
    else:
      return True
    #} end if
  #} end def

  def CheckStatus(self): #{
    fail_msg = \
      "cannot open output file for job number %s" % self.num
    output_file = FileBoxCls(self.output_path, "r", fail_msg)
    self.status = "in progress"
    for output_line in output_file: #{
      if (P2G_SUCCESS == output_line): #{
        self.status = "complete"
        break
      elif ("" == output_line or P2G_FAIL == output_line): #{
        self.status = "failed"
        break
      #} end if
    #} end for
    output_file.close()
    return
  #} end def

  def GetJobLine(self): #{
    fail_msg = "cannot read job file"
    job_file = FileBoxCls(self.path, "r", fail_msg)
    job_line = None
    # the job line should be the last line in the file
    for line in job_file: #{
      job_line = line
    #} end for
    if (None == job_line or "pair_to_genome/calculate" not in line): #{
      raise P2GStatusCheckerError \
        ("cannot get job line from file: %s" % self.path)
    #} end if
    job_file.Close()
    return job_line
  #} end def

  def JobName(self): #{
    #return "%s-p2g" % self.lib
    return os.path.basename(self.path)
  #} end def

  def NewOutPath(self): #{
    if (None == self.output_path): #{
      old_out_num = 0
    else:
      old_out_match = re.search(r"o(?P<out_num>\d+)$", self.output_path)
      if (None == old_out_match): #{
        raise P2GStatusCheckerError("error getting new stdout path for job, "
          "cannot get number from \"%s\"" % self.output_path)
      #} end if
      old_out_num = int(old_out_match.group('out_num'))
    #} end if
    new_out_num = old_out_num + 1
    new_out_path = "%s.o%i" % (self.path, new_out_num)
    CheckNewFilePath(new_out_path, "new stdout")
    return new_out_path
  #} end def

  def NewErrPath(self): #{
    if (None == self.output_path): #{
      old_out_num = 0
    else:
      old_out_match = re.search(r"o(?P<out_num>\d+)$", self.output_path)
      if (None == old_out_match): #{
        raise P2GStatusCheckerError("error getting new stdout path for job, "
          "cannot get number from \"%s\"" % self.output_path)
      #} end if
      old_out_num = int(old_out_match.group('out_num'))
    #} end if
    new_err_num = old_out_num + 1
    new_err_path = "%s.e%i" % (self.path, new_err_num)
    CheckNewFilePath(new_err_path, "new stderr")
    return new_err_path
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class P2GStatusCheckerError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Checks status of pair-to-genome support "
    "cluster jobs.")
  args = [ "JOBS_DIR", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  resubmit_group = OptionGroup(parser, "Resubmission Options")
  resubmit_group.add_option("-r", "--resubmit",
                    action="store_true",
                    help="Resubmit all incomplete jobs.")
  resubmit_group.add_option("-c", "--cluster-head",
                    dest="cluster_head",
                    help="Specify cluster head node if resubmitting jobs to "
                         "cluster")
  resubmit_group.add_option("--mem", "--memory",
                    help="The memory requirement of the resubmitted jobs "
                         "[default:\"%default\"].")
  resubmit_group.add_option("--hostname",
                    help="The hostname(s) to resubmit to "
                         "[default: read from barnacle.cfg].")
  resubmit_group.add_option("--queue",
                    help="The queue(s) to resubmit to [default: read "
                      "from barnacle.cfg].")
  resubmit_group.add_option("-w", "--wall-time",
                    metavar="H:MM:SS",
                    help="When resubmitting, use this value as the "
                         "wall-time option.")
  resubmit_group.add_option("--email",
                    help="E-mail status updates on resubmitted jobs to the "
                         "given email address")
  parser.add_option_group(resubmit_group)
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  misc_group.add_option("--log-file",
                    dest="log_file_name", metavar="FILE",
                    help="Log all messages in FILE")
  misc_group.add_option("-t", "--terse",
                    action="store_true",
                    help="Only write number complete and total number of jobs")
  misc_group.add_option("-q", "--quiet",
                    action="store_true",
                    help="Only write output to log-file, not to the screen.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  misc_group.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(resubmit=False,
                      mem="1G",
                      hostname="q*",
                      queue="all.q",
                      dpt=False,
                      terse=False,
                      quiet=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  path_errors = list()
  CheckDirPath(options.jobs_dir, "pair-to-genome cluster jobs", path_errors)
  if (options.quiet and None == options.log_file_name): #{
    path_errors.append \
      ("please specify a log file to use when using the quiet option")
  #} end if
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
  if (parser.num_args == len(args)): #{
    options.jobs_dir = EnsureAbsPath(args[0])
    if (CheckPaths(options)): #{
      try:
        p2g_status_checker = P2GStatusCheckerCls(options)
        WriteCommand(p2g_status_checker, sys.argv)
        p2g_status_checker.CheckP2GStatus()
        p2g_status_checker.OutputP2GStatus()
        if (options.resubmit): #{
          if (0 < len(p2g_status_checker.unfinished_jobs)): #{
            p2g_status_checker.ResubmitJobs()
          else:
            LogMsg(p2g_status_checker,
              "All jobs are complete: nothing to resubmit")
          #} end def
        #} end if
      except (MyError), e:
        ErrMsg("ERROR while checking pair-to-genome support status:\n"
               "  %s" % e)
        return ES_RUN_ERR
      # end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a path to a directory containing "
      "pair-to-genome cluster jobs for that library (JOBS_DIR).")
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
