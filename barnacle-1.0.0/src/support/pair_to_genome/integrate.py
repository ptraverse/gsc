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
from check_status  import P2GStatusCheckerCls
from group import P2GGroupCls
#from process_candidates.filter import FilterCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "P2G INTEGRATE SUCCESS"
MSG_FAIL = "P2G INTEGRATE FAIL"

class P2GIntegratorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    self.group_parser = None
    self.output_file  = None
    self.num_groups   = 0
    self.groups_without_reads = list()
  #} end def

  def __del__(self): #{
    if (hasattr(self, "output_file") and None != self.output_file): #{
      self.output_file.Close()
    #} end if
    CloseLogFile(self)
  #} end def

  def IntegrateP2GSupport(self): #{
    start_time = time.time()
    LogMsg(self, "Integrating pair-to-genome support for %s..." %
      self.options.lib)
    # check that all pair-to-genome support jobs are complete
    self.CheckP2GStatus()
    # open the input and output files
    self.OpenFiles()
    # loop through the pair-to-genome results files
    for p2g_path in self.GetPairToGenomePaths(): #{
      self.IntegrateP2GFile(p2g_path)
    #} end for
    # output warnings about groups without reads
    self.WarnAboutReadlessGroups()
    LogMsg(self, "Number of groups: %i" % self.num_groups)
    LogMsg(self, "Time spent integrating pair-to-genome support "
      "results: %s" % TimeSpent(start_time))
  #} end def

  def CheckP2GStatus(self): #{
    # setup the status checker options
    status_checker_options = self.options
    status_checker_options.debug = False
    status_checker_options.jobs_dir = self.options.jobs_dir
    status_checker_options.terse = False
    status_checker_options.quiet = False
    # check the status of the pair-to-genome jobs
    p2g_status_checker = \
      P2GStatusCheckerCls(status_checker_options, self.log_info)
    p2g_status_checker.CheckP2GStatus()
    if (self.log_info['debug']): #{
      p2g_status_checker.OutputP2GStatus()
    #} end if
    # ensure that the status is "complete"
    if ("complete" == p2g_status_checker.status): #{
      LogMsg(self, "All pair-to-genome support jobs complete")
    elif ("failed" == p2g_status_checker.status):
      raise P2GIntegratorError \
        ("pair-to-genome support jobs failed")
    elif ("in progress" == p2g_status_checker.status):
      raise P2GIntegratorError \
        ("pair-to-genome support jobs incomplete")
    elif ("unknown"  == p2g_status_checker.status):
      raise P2GIntegratorError \
        ("cannot determine status of pair-to-genome support jobs: %s" %
         p2g_status_checker.status)
    #} end if
  #} end if

  def OpenFiles(self): #{
    # create a group parser
    self.group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    # create the output file
    fail_msg = "could not create output file"
    self.output_file = FileBoxCls(self.options.out_path, "w", fail_msg)
  #} end def

  def GetPairToGenomePaths(self): #{
    LogMsg(self, "Getting pair-to-genome result file paths...")
    p2g_paths = list()
    job_dir_template = os.path.join(self.options.jobs_dir, "job_*")
    job_dirs = glob(job_dir_template)
    for job_dir in sorted(job_dirs, key=GetJobNum): #{
      job_num = GetJobNum(job_dir)
      p2g_path = os.path.join(job_dir, "%s.%i.p2g.out" %
        (self.options.lib, job_num))
      CheckFilePath(p2g_path, "pair-to-genome results")
      p2g_paths.append(p2g_path)
      if (len(p2g_paths) != job_num): #{
        raise P2GIntegratorError("pair-to-genome results file missing: %i" %
          len(p2g_paths))
      #} end if
    #} end for
    LogMsg(self, "Integrating pair-to-genome result files...")
    return p2g_paths
  #} end def

  def IntegrateP2GFile(self, p2g_path): #{
    DebugMsg(self, "Integrating pair-to-genome file: %s" % p2g_path)
    group = None
    fail_msg = "cannot open pair-to-genome results file"
    p2g_file = FileBoxCls(p2g_path, "r", fail_msg)
    for p2g_line in p2g_file: #{
      DebugMsg(self, "LINE: %s" % p2g_line)
      # count the group
      self.num_groups += 1
      # parse the pair-to-genome line
      p2g_support = P2GGroupCls(self.options, self.log_info)
      p2g_support.ParseSupportString(p2g_line)
      # check that the group had some reads at least
      if (1 > p2g_support.num_reads): #{
        self.groups_without_reads.append("%i" % p2g_support.group_id)
      #} end if
      # get a group from the groups file
      if (None == group or
          p2g_support.group_id > group.id):
        try: #{
          DebugMsg(self, "Getting next group...")
          group = self.group_parser.GetNextGroup()
        except StopIteration:
          raise P2GIntegratorError \
            ("Unexpected end of groups file: %s\n  while integrating: %s" %
             (self.group_parser.data_file_path, p2g_path))
        #} end try
      # allow for groups having been removed from the groups file
      if (p2g_support.group_id < group.id): #{
        continue
      #} end if
      # ensure that the group ids match up
      if (p2g_support.group_id != group.id): #{
        raise P2GIntegratorError("Inconsistent group ids: %i from %s, " %
          (p2g_support.group_id, p2g_path) +
          "%i from %s" % (group.id, self.options.barnacle_path))
      #} end if
      # add the pair-to-genome support to the group
      self.AddSupportToGroup(group, p2g_support)
      # apply any pair-to-genome filters given
      #self.ApplyFilters(group)
      # write the group to the new output file(s)
      self.WriteGroup(group)
    #} end for
    p2g_file.close()
  #} end def

  def AddSupportToGroup(self, group, p2g_support): #{
    DebugMsg(self, "Adding support to group %i" % group.id)
    p2g_errors = list()
    if (p2g_support.p2g_exonic > p2g_support.p2g_all): #{
      p2g_errors.append("exonic: %i cannot be more than total: %i" %
        (p2g_support.p2g_exonic, p2g_support.p2g_all))
    #} end if
    if (p2g_support.p2g_filt_exonic > p2g_support.p2g_filt_all): #{
      p2g_errors.append("filtered exonic: %i " % p2g_support.p2g_filt_exonic +
        "cannot be more than filtered total: %i" % p2g_support.p2g_filt_all)
    #} end if
    if (p2g_support.p2g_filt_exonic > p2g_support.p2g_exonic): #{
      p2g_errors.append("filtered exonic: %i " % p2g_support.p2g_filt_exonic +
        "cannot be more than unfiltered: %i" % p2g_support.p2g_exonic)
    #} end if
    if (p2g_support.p2g_filt_all > p2g_support.p2g_all): #{
      p2g_errors.append("filtered total: %i " % p2g_support.p2g_filt_all +
        "cannot be more than unfiltered: %i" % p2g_support.p2g_all)
    #} end if
    if (0 < len(p2g_errors)): #{
      raise P2GIntegratorError("invalid pair-to-genome support values:\n"
        "  %s" % "\n  ".join(p2g_errors))
    #} end if
    group.pair_to_genome_exonic          = p2g_support.p2g_exonic
    group.pair_to_genome_all             = p2g_support.p2g_all
    group.pair_to_genome_filtered_exonic = p2g_support.p2g_filt_exonic
    group.pair_to_genome_filtered_all    = p2g_support.p2g_filt_all
  #} end def

  def ApplyFilters(self, group): #{
    # setup the filtering options
    filter_options = self.options
    filter_options.sort_by_r2cu = False
    filter_options.gene_names_path = None
    # create the filtering object
    filter = QuickFilterCls(filter_options, self.log_info)
    # apply the filters
    filter.ApplyPairToGenomeFilter(group)
    DebugMsg(self, "Group %i: p2g=%i %s" %
      (group.id, group.pair_to_genome_exonic, group.Status()))
  #} end def

  def WriteGroup(self, group): #{
    DebugMsg(self, "Writing group %i" % group.id)
    # write the group to the all groups file
    self.output_file.Write(group.FullDataString())
    # if the group passed the filters
    #if (group.PassedFilters()): #{
    #  # count the passing group
    #  self.num_passing += 1
    #  # write the group to the all passing file
    #  self.output_files['passing'].Write(group_string)
    #  # write the group to the appropriate type file
    #  #self.output_files[group.Type()].Write(group_string)
    #} end if
  #} end def

  def WarnAboutReadlessGroups(self): #{
    if (0 == len(self.groups_without_reads)): #{
      return
    #} end if
    LogMsg(self, "WARNING: SAMtools view returned no reads for %i groups" %
      len(self.groups_without_reads))
    DebugMsg(self, "  Groups without reads: %s" %
      ",".join(self.groups_without_reads))
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class P2GIntegratorError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Merges split pair-to-genome support results "
    "and integrates them into the Barnacle results")
  args = [ "LIB", "BARNACLE_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(
                      dpt=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  path_errors = list()
  # check the groups path
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  # get the input directory
  (input_dir, file_name) = os.path.split(options.barnacle_path)
  # get and check the jobs directory
  options.jobs_dir = os.path.join(os.path.dirname(input_dir),
    "cluster_p2g")
  CheckDirPath(options.jobs_dir, "pair-to-genome jobs",
    path_errors, create=False)
  # get and check the output path
  options.output_dir = GetOutDir(input_dir, "with_p2g")
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors, create=True)
    # check that the output file does not already exist
    options.out_path = os.path.join(options.output_dir, file_name)
    CheckNewFilePath(options.out_path, "pair-to-genome support", path_errors)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "p2g_integrate", options.output_dir)
  #} end if
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
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
    options.lib           = args[0]
    options.barnacle_path = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try: #{
        p2g_integrator = P2GIntegratorCls(options)
        WriteCommand(p2g_integrator, sys.argv)
        p2g_integrator.IntegrateP2GSupport()
      except (P2GIntegratorError), e:
        ErrMsg("ERROR while integrating pair-to-genome support:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
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
