#! /usr/bin/env python
"""
integrate.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
from glob import glob
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
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.submit import GetJobNum
from utils.general import (SetupMainClass, TimeSpent, WriteCommand)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, CheckNewFilePath, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls
from check_status  import R2CStatusCheckerCls
from group import R2CGroupCls, R2CMemberCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "R2C INTEGRATE SUCCESS"
MSG_FAIL = "R2C INTEGRATE FAIL"

class R2CIntegratorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    self.group_parser = None
    # group_to_jobs_map[group_id] = list of job indices
    self.group_to_jobs_map = dict()
    self.r2c_files    = list()
    self.output_file  = None
  #} end def

  def __del__(self): #{
    if (hasattr(self, "output_file") and None != self.output_file): #{
      self.output_file.Close()
    #} end if
    CloseLogFile(self)
  #} end def

  def IntegrateR2CSupport(self): #{
    LogMsg(self, "Integrating read-to-contig support for %s..." %
      self.options.lib)
    start = time.time()
    # check that all read-to-contig support jobs are complete
    self.CheckR2CStatus()
    # open the input and output files
    self.Setup()
    # loop through the grouped candidates file
    for group in self.group_parser: #{
      DebugMsg(self, "Integrating read-to-contig support for group %i" %
        group.id)
      self.IntegrateR2CForGroup(group)
      self.WriteGroup(group)
    #} end for
    self.output_file.Close()
    LogMsg(self, "Time spent integrating read-to-contig support "
      "results: %s" % TimeSpent(start))
  #} end def

  def CheckR2CStatus(self): #{
    # setup the status checker options
    status_checker_options = self.options
    status_checker_options.debug = False
    status_checker_options.jobs_dir = self.options.jobs_dir
    status_checker_options.terse = False
    status_checker_options.quiet = False
    status_checker_options.hostname = None
    status_checker_options.queue = None
    # check the status of the read-to-contig jobs
    r2c_status_checker = R2CStatusCheckerCls(status_checker_options,
        self.log_info)
    r2c_status_checker.CheckR2CStatus()
    if (self.log_info['debug']): #{
      r2c_status_checker.OutputR2CStatus()
    #} end if
    # ensure that the status is "complete"
    if ("complete" == r2c_status_checker.status): #{
      LogMsg(self, "All read-to-contig support jobs complete")
    elif ("failed" == r2c_status_checker.status):
      raise R2CIntegratorError \
        ("read-to-contig support jobs failed")
    elif ("in progress" == r2c_status_checker.status):
      raise R2CIntegratorError \
        ("read-to-contig support jobs incomplete")
    elif ("unknown"  == r2c_status_checker.status):
      raise R2CIntegratorError \
        ("cannot determine status of read-to-contig support jobs: %s" %
         r2c_status_checker.status)
    #} end if
  #} end if

  def Setup(self): #{
    # create a group parser
    self.group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    # read in the group-to-job map
    DebugMsg(self, "Loading group-to-jobs map...")
    fail_msg = "cannot open group-to-jobs map file"
    group_to_jobs_file = FileBoxCls(self.options.group_to_jobs_path,
      "r", fail_msg)
    for group_to_jobs_line in group_to_jobs_file: #{
      (group_id_str, job_indices_str) = group_to_jobs_line.split(":")
      group_id = int(group_id_str)
      if (group_id in self.group_to_jobs_map): #{
        raise R2CIntegratorError("duplicated group id in group-to-jobs "
          "map: %i" % group_id)
      #} end if
      if ("" != job_indices_str): #{
        self.group_to_jobs_map[group_id] = map(int, job_indices_str.split(","))
      #} end if
      #DebugMsg(self, "Mapping group %i to jobs %s" % (group_id,
      #  ",".join(map(str, self.group_to_jobs_map[group_id]))))
    #} end for
    group_to_jobs_file.Close()
    # open the results files
    DebugMsg(self, "Opening read-to-contig support results files")
    job_dir_template = os.path.join(self.options.jobs_dir, "job_*")
    job_dirs = glob(job_dir_template)
    for job_dir in sorted(job_dirs, key=GetJobNum): #{
      job_num = GetJobNum(job_dir)
      r2c_path = os.path.join(job_dir, "%s.%i.r2c.out" %
        (self.options.lib, job_num))
      CheckFilePath(r2c_path, "pair-to-genome results")
      r2c_file = R2CResultsFileCls(r2c_path, log_info=self.log_info)
      self.r2c_files.append(r2c_file)
      if (len(self.r2c_files) != job_num): #{
        raise P2GIntegratorError("pair-to-genome results file missing: %i" %
          len(self.r2c_files))
      #} end if
    #} end for
    # create the output file
    fail_msg = "could not create output file"
    self.output_file = FileBoxCls(self.options.out_path, "w", fail_msg)
  #} end def

  def IntegrateR2CForGroup(self, group): #{
    if (group.id not in self.group_to_jobs_map): #{
      DebugMsg(self, "Skipping group %i: not in group_to_jobs_map: " %
        group.id)
      return
    #} end if
    support = R2CGroupCls(group.id, log_info=self.log_info)
    for job_index in self.group_to_jobs_map[group.id]: #{
      #DebugMsg(self, "Group: %i, Job: %i" % (group.id, job_index))
      self.AddSupportFromFile(support, job_index)
    #} end for
    for member in group.members: #{
      self.CalculateSupport(support, member)
    #} end for
  #} end def

  def AddSupportFromFile(self, curr_group, job_index): #{
    # external directories use base 1, convert to base 0
    curr_file = self.r2c_files[job_index-1]
    DebugMsg(self, "Adding support to group %i from file %s" %
      (curr_group.group_id, curr_file.file.path))
    # skip groups in the current results file that
    # come before the current group
    while (curr_file.BeforeGroup(curr_group.group_id)): #{
      if (None != curr_file.curr_member): #{
        DebugMsg(self, "Skipping support for member %s" %
          curr_file.curr_member.member_id)
      #} end if
      curr_file.GetMember()
    #} end while
    if (curr_file.integrated): #{
      DebugMsg(self, "File has been fully integrated")
      return
    #} end if
    DebugMsg(self, "Target group: %i, Current member group: %i" %
      (curr_group.group_id, curr_file.curr_member.group_id))
    # add members for the current group to the current group
    while (curr_file.GroupIsCurrent(curr_group.group_id)): #{
      DebugMsg(self, "Adding support for member %s" %
        curr_file.curr_member.member_id)
      # if the group already has the current member,
      # add the new support values
      curr_group.AddSupport(curr_file.curr_member)
      # get the next member
      curr_file.GetMember()
    #} end while
  #} end def

  def CalculateSupport(self, support, member): #{
    DebugMsg(self, "Calculating support for group member %s, (%s)" %
      (member.IDString(), ",".join(support.members.keys())))
    if (member.IDString() in support.members): #{
      member_support = support.members[member.IDString()]
      member_support.CalculateSupport()
      member.read_to_contig_all     = member_support.min_support
      member.avg_read_to_contig_all = member_support.average_support
      member.read_to_ctg_unique     = member_support.min_support_unique
      member.avg_read_to_ctg_unique = member_support.average_support_unique
    else:
      member.read_to_contig_all     = 0
      member.avg_read_to_contig_all = 0
      member.read_to_ctg_unique     = 0
      member.avg_read_to_ctg_unique = 0
    #} end if
  #} end def

  def WriteGroup(self, group): #{
    DebugMsg(self, "Writing group %i" % group.id)
    # write the group to the all groups file
    self.output_file.Write(group.FullDataString())
  #} end def
#} end class

class R2CResultsFileCls: #{
  def __init__(self, path, log_info=None): #{
    self.log_info = log_info
    fail_msg = "cannot open read-to-contig support results file"
    self.file = FileBoxCls(path, "r", fail_msg)
    self.integrated = False
    self.curr_member = None
  #} end def

  def __del__(self): #{
    self.file.Close()
  #} end def

  def BeforeGroup(self, group_id): #{
    if (self.integrated): #{
      return False
    #} end if
    if (None == self.curr_member or self.curr_member.group_id < group_id): #{
      return True
    #} end if
    return False
  #} end def

  def GroupIsCurrent(self, group_id): #{
    if (self.integrated or None == self.curr_member): #{
      return False
    #} end if
    if (self.curr_member.group_id == group_id): #{
      return True
    #} end if
    return False
  #} end def

  def GetMember(self): #{
    if (self.integrated): #{
      DebugMsg(self, "Not getting member, file already fully integrated.")
      return
    #} end if
    DebugMsg(self, "Getting member...")
    try: #{
      member_line = self.file.next()
      # create a new member from the current line, store it as "curr_member"
      (member_id, support_list) = member_line.split(" ")
      self.curr_member = R2CMemberCls(member_id, log_info=self.log_info)
      # store support values
      self.curr_member.InitializeSupport(support_list)
      DebugMsg(self, "New member: %s" % self.curr_member.DebugString())
    except StopIteration:
      DebugMsg(self, "Integrated all support from %s" % self.file.path)
      self.curr_member = None
      self.integrated = True
      self.file.Close()
      return
    #} end try
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class R2CIntegratorError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Merges split read-to-contig support results file "
    "and integrates it into the chimeric transcript results")
  args = [ "LIB", "BARNACLE_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--r2c-dir",
                    dest = "jobs_dir", metavar="DIR",
                    help="Use this directory for the read-to-contig support "
                         "results, rather than the default.")
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
  CheckFilePath(options.barnacle_path, "groups", path_errors)
  # get the input directory
  (input_dir, file_name) = os.path.split(options.barnacle_path)
  # get and check the jobs directory
  if (None == options.jobs_dir): #{
    options.jobs_dir = os.path.join(os.path.dirname(input_dir),
      "cluster_r2c")
  #} end if
  CheckDirPath(options.jobs_dir, "read-to-contig jobs",
    path_errors, create=False)
  # get and check the group-to-jobs map path
  options.group_to_jobs_path = os.path.join(options.jobs_dir,
    "group2jobs_map")
  CheckFilePath(options.group_to_jobs_path, "group-to-jobs map", path_errors)
  # get and check the output path
  options.output_dir = GetOutDir(input_dir, "with_r2c")
  CheckDirPath(options.output_dir, "output", path_errors, create=True)
  # check that the output file does not already exist
  options.out_path = os.path.join(options.output_dir, file_name)
  CheckNewFilePath(options.out_path, "read-to-contig support", path_errors)
  # get the log file name
  options.log_file_name = GetLogPath(options.barnacle_path,
    "r2c_integrate", options.output_dir)
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
      try: #{
        r2c_integrator = R2CIntegratorCls(options)
        WriteCommand(r2c_integrator, sys.argv)
        r2c_integrator.IntegrateR2CSupport()
      except (MyError), e:
        ErrMsg("ERROR while integrating read-to-contig support:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); and the path to a "
      "BARNACLE group data file for that library (BARNACLE_FILE).")
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
