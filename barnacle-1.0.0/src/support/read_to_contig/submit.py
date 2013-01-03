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
from utils import submit
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand, GetOptions,
  CheckConfigCommands, GetCommand, GetClusterValue)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, CheckNewFilePath, EnsureDirectoryExists, FileBoxCls)
from utils.subprocesses import (RunCommandFromString, SortFile)
from parsers.candidate_group_parser import CandidateGroupParserCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "R2C SUBMIT SUCCESS"
MSG_FAIL = "R2C SUBMIT FAIL"

class R2CSubmitterCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    if (None == self.options.cluster_head): #{
      local = True
    else:
      local = False
    #} end if
    CheckConfigCommands(self, "python", local)
    self.jobs = R2CJobCollectionCls(options, log_info=self.log_info)
    # group_to_jobs_map is a file
    self.group_to_jobs_map = None
    self.short_ctgs = set()
    if (None == self.options.hostname): #{
      self.options.hostname = GetClusterValue(self, "hostname")
    #} end if
    if (None == self.options.queue): #{
      self.options.queue = GetClusterValue(self, "queue")
    #} end if
  #} end def

  def __del__(self): #{
    if (hasattr(self, "group_to_jobs_map") and
        None != self.group_to_jobs_map): #{
      self.group_to_jobs_map.Close()
    #} end if
    for job in self.jobs: #{
      job.Close()
    #} end for
    CloseLogFile(self)
  #} end def

  def Run(self): #{
    LogMsg(self, "Creating parallel read-to-contig jobs...")
    start = time.time()
    self.Setup()
    LogMsg(self, "Processing groups...")
    process_start = time.time()
    for group in self.groups_file: #{
      DebugMsg(self, "Processing group %i..." % group.id)
      DebugMsg(self, group.FullDataString().rstrip())
      group_start = time.time()
      try: #{
        self.ProcessGroup(group)
      except MyError, e:
        raise R2CSubmitterError("Error processing group %i:\n%s" %
          (group.id, e))
      #} end try
      DebugMsg(self, "Time spent processing group: %s\n" %
        TimeSpent(group_start) + "--------------------------\n")
    #} end for
    LogMsg(self, "Time spent processing groups: %s" %
      TimeSpent(process_start))
    if (None == self.options.cluster_head): #{
      self.RunJobsLocally()
    else:
      jobs_path = self.CreateJobsFile()
      self.SubmitJobs(jobs_path)
    #} end if
    if (self.options.no_short_ctgs and 0 < len(self.short_ctgs)): #{
      LogMsg(self, "WARNING: %i contigs shorter than the read-length have "
        "been ignored." % len(self.short_ctgs))
    #} end if
    LogMsg(self, "Time spent creating parallel read-to-contig "
      "jobs: %s" % TimeSpent(start))
  #} end def

  def Setup(self): #{
    self.groups_file = CandidateGroupParserCls(self.options.barnacle_path)
    fail_msg = "cannot create groups-to-jobs map file"
    groups_to_jobs_path = os.path.join(self.options.cluster_dir,
      "group2jobs_map")
    self.group_to_jobs_map = FileBoxCls(groups_to_jobs_path, "w", fail_msg)
  #} end def

  def ProcessGroup(self, group): #{
    job_indices = set()
    member_strings = list()
    # get member strings for group
    for member in group.members: #{
      member_strings.append(self.GetMemberString(member))
    #} end for
    # write member strings to appropriate jobs files
    for member_index in range(len(member_strings)): #{
      mate_string_list = (member_strings[:member_index] +
        member_strings[member_index+1:])
      if (0 == len(mate_string_list)): #{
        mates_string = "None"
      else:
        mates_string = ";".join(mate_string_list).replace(" ", ",")
      #} end if
      full_member_str = member_strings[member_index] + \
        " MATES:%s" % mates_string
      DebugMsg(self, "  MEMBER STRING: %s" % full_member_str)
      ctg_id = group.members[member_index].contig_info.id
      ctg_len = group.members[member_index].contig_info.length
      if (self.options.no_short_ctgs and
          self.options.read_length > ctg_len): #{
        self.short_ctgs.add(ctg_id)
        DebugMsg(self, "Not adding member to job: contig shorter "
          "than read-length.")
        continue
      #} end if
      job_index = self.jobs.ProcessMember(ctg_id, full_member_str)
      job_indices.add(str(job_index))
    #} end for
    jobs_string = "%i:%s" % (group.id, ",".join(job_indices))
    self.group_to_jobs_map.WriteLine(jobs_string)
  #} end def

  def GetMemberString(self, member): #{
    ExtremeDebugMsg(self, "  Getting Member String for %s..." %
      member.IDString())
    (left, right) = self.GetBufferedEventCoords(member)
    member_info = [
      ConvertContigID(member.contig_info.id),
      member.IDString(),
      "%i-%i" % (left, right),
    ]
    return " ".join(member_info)
  #} end def

  def GetBufferedEventCoords(self, member): #{
    #if (member.gap): #{
    #  # if the gap is internal to the contig
    #  if (member.GapIsInternal()): #{
    #    left  = member.align_info_B.ctg_start
    #    right = member.align_info_B.ctg_end
    #  else:
    #    # figure out which coordinates to use
    #    (realign_start, realign_end,
    #     dup_start, dup_end) = member.RealignmentCoords()
    #    #realign_start = member.align_info_B.ctg_start
    #    #realign_end   = member.align_info_B.ctg_end
    #    #if ('dup_start' in member.meta_fields): #{
    #    #  (olapA, fractionA) = CalcOverlap(member.meta_fields['gap_start'],
    #    #    member.meta_fields['gap_end'], realign_start, realign_end)
    #    #  (olapB, fractionB) = CalcOverlap(member.meta_fields['gap_start'],
    #    #    member.meta_fields['gap_end'], member.meta_fields['dup_start'],
    #    #    member.meta_fields['dup_end'])
    #    #  ctg_str = "ctg:%i-%i(%i)" % (member.align_info_B.ctg_start,
    #    #    member.align_info_B.ctg_end, olapA)
    #    #  dup_str = "dup:%i-%i(%i)" % (member.meta_fields['dup_start'],
    #    #    member.meta_fields['dup_end'], olapB)
    #    #  DebugMsg(self, "Checking which coordinates to use for edge-gap\n"
    #    #    "gap:%i-%i %s %s" % (member.meta_fields['gap_start'],
    #    #     member.meta_fields['gap_end'], ctg_str, dup_str))
    #    #  if (olapB > olapA): #{
    #    #    realign_start = member.meta_fields['dup_start']
    #    #    realign_end   = member.meta_fields['dup_end']
    #    #    DebugMsg(self, "Using DUP coordinates")
    #    #  else:
    #    #    DebugMsg(self, "Using CTG coordinates")
    #    #  #} end if
    #    #} end if
    #    # if the gap is at the start of the contig
    #    if (member.StartOfContig() == member.meta_fields['gap_start']): #{
    #      left  = min(realign_end, member.meta_fields['gap_end'])
    #      right = max(realign_end, member.meta_fields['gap_end'])
    #    # if the gap is at the end of the contig
    #    else:
    #      left  = min(realign_start, member.meta_fields['gap_start'])
    #      right = max(realign_start, member.meta_fields['gap_start'])
    #    #} end if
    #  #} end if
    #else:
    #  left  = min(member.align_info_A.ctg_end, member.align_info_B.ctg_start)
    #  right = max(member.align_info_A.ctg_end, member.align_info_B.ctg_start)
    #} end if
    # get the standard event region
    region = member.ContigEventRegion()
    ExtremeDebugMsg(self, "    Before expansion: %s" % region)
    # widen the region to satisfy the minimum overlap requirement and
    # ensure that the group does not extend past the edges of the contig
    min_left  = 2
    max_right = member.contig_info.length-1
    #buffered_left  = max(left  - self.options.min_overlap, min_left)
    #buffered_right = min(right + self.options.min_overlap, max_right)
    region.min = max(region.min - self.options.min_overlap, min_left)
    region.max = min(region.max + self.options.min_overlap, max_right)
    #if ((buffered_right - buffered_left) > 1000): #{}
    if (1000 < region.Span()): #{
      if (member.gap): #{
        if (member.GapIsInternal()): #{
          event_type = "internal gap"
        else:
          event_type = "edge gap"
        #} end if
      else:
        event_type = "split"
      #} end if
      LogMsg(self, "WARNING: %s has wide %s event coordinates: %i" %
        (member.IDString(), event_type, region.Span()))
    #} end if
    # TODO: split wide gap regions into two regions?
    return (region.min, region.max)
  #} end def

  def CreateJobsFile(self): #{
    jobs_path = os.path.join(self.options.cluster_dir, "read2ctg_jobs")
    LogMsg(self, "Creating jobs file: %s" % jobs_path)
    fail_msg = "cannot create jobs file"
    jobs_file = FileBoxCls(jobs_path, "w", fail_msg)
    sort_time_total = 0
    for job in self.jobs: #{
      job.Close()
      sort_start = time.time()
      job.SortFileByContig()
      DebugMsg(self, "Time spent sorting job data: %s" % TimeSpent(sort_start))
      sort_time_total += (time.time() - sort_start)
      jobs_file.WriteLine(job.Command(local=False))
    #} end for
    sort_time_avg = sort_time_total / len(self.jobs)
    LogMsg(self, "Mean time spent sorting job data: %s" %
      TimeSpent(sort_time_avg, calc_elapsed=False))
    jobs_file.Close()
    return jobs_file.path
  #} end def

  def RunJobsLocally(self): #{
    LogMsg(self, "Running jobs locally...")
    for job in self.jobs: #{
      job.Close()
      sort_start = time.time()
      job.SortFileByContig()
      LogMsg(self, "Time spent sorting job data: %s" % TimeSpent(sort_start))
      fail_msg = "cannot create job command file"
      # create job-file in case cluster resubmission occurs
      job_file = FileBoxCls(job.CommandPath(), "w", fail_msg)
      job_file.WriteLine(job.Command(local=False))
      job_file.Close()
      fail_msg = "cannot create stdout file for job"
      out_file = {'path':job.StdOutPath(), 'mode':"w", 'fail_msg':fail_msg}
      fail_msg = "cannot create stderr file for job"
      err_file = {'path':job.StdErrPath(), 'mode':"w", 'fail_msg':fail_msg}
      job_command = job.Command(local=True)
      DebugMsg(self, job_command)
      return_code = RunCommandFromString(job_command,
        stdout=out_file, stderr=err_file, dpt=self.options.dpt)
      if (0 != return_code and not self.log_info['debug']): #{
        LogMsg(self, job_command)
      #} end if
      if (0 > return_code): #{
        raise R2CSubmitterError("Read-to-contig job command was terminated "
          "by signal %i" % return_code)
      elif (0 < return_code):
        raise R2CSubmitterError("Error running read-to-contig job "
          "command: %i" % return_code)
      #} end if
    #} end for
  #} end def

  def SubmitJobs(self, jobs_path): #{
    LogMsg(self, "Submitting read-to-contig jobs...")
    submit_options = GetOptions(submit)
    submit.SetSubmitOptions(submit_options, "%s-r2c" % self.options.lib,
      jobs_path, False, self.options)
    if (submit.CheckPaths(submit_options)): #{
      submitter = submit.JobSubmitterCls(submit_options,
        log_info=self.log_info)
      submitter.Submit()
    else:
      raise R2CSubmitterError("cannot submit jobs: "
        "errors in submit options")
    #} end if
  #} end def
#} end class

class R2CJobCollectionCls: #{
  def __init__(self, options, log_info=None): #{
    self.options  = options
    self.log_info = log_info
    # internally job list uses base 0, but external directories use base 1
    # job_lookup[ctg_id] = job_index
    self.job_lookup = dict()
    self.jobs = list()
  #} end def

  def __iter__(self): #{
    return self.jobs.__iter__()
  #} end def

  def __len__(self): #{
    return len(self.jobs)
  #} end def

  def ProcessMember(self, ctg_id, member_str): #{
    DebugMsg(self, "Processing member with contig: %s" % ctg_id)
    if (ctg_id not in self.job_lookup): #{
      self.AddNewContig(ctg_id)
    #} end if
    job_index = self.job_lookup[ctg_id]
    self.jobs[job_index].ProcessMember(member_str)
    # internally job list uses base 0, but external directories use base 1
    return self.jobs[job_index].index
  #} end def

  def AddNewContig(self, ctg_id): #{
    DebugMsg(self, "Adding new contig to job collection")
    if (0 == len(self.jobs) or self.jobs[-1].IsFull()): #{
      DebugMsg(self, "Adding new job to job collection, num jobs: %i" %
        (len(self.jobs)+1))
      new_job = R2CJobCls(len(self.jobs)+1, self.options,
        log_info=self.log_info)
      self.jobs.append(new_job)
    #} end if
    self.jobs[-1].AddNewContig(ctg_id)
    # internally job list uses base 0, but external directories use base 1
    self.job_lookup[ctg_id] = self.jobs[-1].index-1
  #} end def
#} end class

class R2CJobCls: #{
  def __init__(self, index, options, log_info=None): #{
    self.options  = options
    self.log_info = log_info
    self.index    = index
    self.contigs  = set()
    self.file     = None
  #} end def

  def IsFull(self): #{
    num_contigs = len(self.contigs)
    if (num_contigs < self.options.ctgs_per_job): #{
      return False
    #} end if
    return True
  #} end def

  def AddNewContig(self, ctg_id): #{
    self.contigs.add(ctg_id)
  #} end def

  def CreateJobDataFile(self): #{
    self.job_dir = os.path.join(self.options.cluster_dir,
      "job_%i" % self.index)
    EnsureDirectoryExists(self.job_dir)
    job_path = os.path.join(self.job_dir, "%s.%i.r2c.data" %
      (self.options.lib, self.index))
    fail_msg = "cannot create job data file"
    self.file = FileBoxCls(job_path, "w", fail_msg)
  #} end def

  def ProcessMember(self, member_str): #{
    if (None == self.file): #{
      self.CreateJobDataFile()
    #} end if
    self.file.WriteLine(member_str)
  #} end def

  def SortFileByContig(self): #{
    if (None == self.file): #{
      raise R2CSubmitterError \
        ("job data file not created for job %i" % self.index)
    #} end if
    self.Close()
    DebugMsg(self, "Sorting job %i data" % self.index)
    SortFile(self, self.file.path, 1, numeric=False)
  #} end def

  def Command(self, local=False): #{
    python_path = GetCommand(self, "python", local)
    script_path = os.path.join(os.environ['BARNACLE_PATH'],
      "support", "read_to_contig", "calculate.py")
    CheckFilePath(script_path, "read-to-contig support script")
    command_list = ["time", python_path, script_path, self.options.lib,
      self.file.path, self.options.r2c_path, "--min-overlap",
      "%i" % self.options.min_overlap]
    if (self.log_info['debug']): #{
      command_list.append("--debug")
    #} end if
    command = " ".join(command_list)
    if (local): #{
      if ("r2c" in os.path.basename(self.file.path)): #{
        log_type = "calculate"
      else:
        log_type = "r2c.calculate"
      #} end if
      log_path = GetLogPath(self.file.path, log_type, self.options.cluster_dir)
      command += " --log-file %s" % log_path
      return command
    #} end if
    full_command_list = [
      "export BARNACLE_PATH=%s" % os.environ["BARNACLE_PATH"],
      "export PATH=$BARNACLE_PATH:$PATH",
      "export PYTHONPATH=.:$BARNACLE_PATH:$PYTHONPATH",
      command,
    ]
    full_command = ";".join(full_command_list)
    return full_command
  #} end def

  def CommandPath(self): #{
    return os.path.join(self.options.cluster_dir,
      "%s-r2c.%i.sh" % (self.options.lib, self.index))
  #} end def

  def StdOutPath(self): #{
    return "%s.o1" % self.CommandPath()
  #} end def

  def StdErrPath(self): #{
    return "%s.e1" % self.CommandPath()
  #} end def

  def Close(self): #{
    if (None != self.file): #{
      self.file.Close()
    #} end if
  #} end def
#} end class

def ConvertContigID(ctg_id): #{
  # replace ":", "+", and "-" with "_", removing commas and trailing "_"
  return re.sub(r"_$", "", re.sub(r"[:+-],?", "_", ctg_id))
#} end def

#### EXCEPTION CLASSES ####
class R2CSubmitterError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Create and submit parallel jobs for calculating "
    "read-to-contig support")
  args = [ "LIB", "BARNACLE_FILE", "READ_TO_CONTIG_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--min-overlap",
                    type="int", metavar="N",
                    help="Require that reads overlap breakpoints by at "
                         "least N bp. [default: %default]")
  parser.add_option("--ctgs-per-job",
                    type="int", metavar="N",
                    help="Process N contigs in each job. [default=%default]")
  parser.add_option("--short-ctgs",
                    action="store_false", dest="no_short_ctgs",
                    help="When creating read-to-conting support jobs, create "
                      "jobs for all contigs, even those shorter than the "
                      "read length.")
  parser.add_option("--no-short-ctgs",
                    action="store_true",
                    help="When creating read-to-conting support jobs, do not "
                      "create jobs for contigs shorter than the read length. "
                      "[default]")
  parser.add_option("--read-length",
                    type="int", metavar="N",
                    help="Reads are N bases long. This option is required "
                      "when using --no-short-ctgs.")
  parser.add_option("-c", "--cluster-head",
                    dest="cluster_head",
                    help="Cluster head node to submit to.")
  parser.add_option("--mem", "--memory",
                    help="The memory requirement of the jobs "
                         "[default:\"%default\"].")
  parser.add_option("--hostname",
                    help="The hostname(s) to submit to "
                         "[default:\"%default\"].")
  parser.add_option("--queue",
                    help="The queue(s) to submit to [default:\"%default\"].")
  parser.add_option("-w", "--wall-time",
                    metavar="H:MM:SS",
                    help="The maximum time to spend on the job.")
  parser.add_option("--email",
                    help="E-mail status updates on submitted jobs to the "
                         "given email address")
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.set_defaults(min_overlap=5,
                      ctgs_per_job=1500,
                      no_short_ctgs=True,
                      read_length=75,
                      mem="5G",
                      dpt=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (options.no_short_ctgs and None == options.read_length): #{
    ErrMsg("You must provide the read-length when using the "
      "--no-short-ctgs option.")
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  (input_dir, input_file_name) = os.path.split(options.barnacle_path)
  options.cluster_dir = os.path.join(os.path.dirname(input_dir),
    "cluster_r2c")
  CheckFilePath(options.r2c_path, "read-to-contigs alignments", path_errors)
  if (0 == len(path_errors)): #{
    # attempt to create output directory, if it does not already exist
    CheckDirPath(options.cluster_dir, "cluster directory", path_errors,
      create=True)
    options.log_file_name = GetLogPath(options.barnacle_path,
      "r2c.submit", options.cluster_dir)
  #} end if
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
    options.lib           = args[0]
    options.barnacle_path = EnsureAbsPath(args[1])
    options.r2c_path      = EnsureAbsPath(args[2])
    if (CheckPaths(options)): #{
      try: #{
        r2c_submitter = R2CSubmitterCls(options)
        WriteCommand(r2c_submitter, sys.argv)
        r2c_submitter.Run()
      except (MyError), e:
        ErrMsg("ERROR while creating parallel read-to-contig jobs:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle data file for that library (BARNACLE_FILE); and the path to "
      "a bam-formatted read to contig alignment file (READ_TO_CONTIG_FILE).")
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
