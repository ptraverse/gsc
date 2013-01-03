#! /usr/bin/env python
"""
submit.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, re, sys, time, traceback

# import custom modules
from version import VERSION
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  CheckConfigCommands)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir)
from parsers.candidate_group_parser import CandidateGroupParserCls
from jobs import HybridMetaJobCollectionCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

class HybridSubmitterCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    if (None == self.options.cluster_head): #{
      local = True
    else:
      local = False
    #} end if
    CheckConfigCommands(self, "python", local)
    self.job_collections = HybridMetaJobCollectionCls(options,
      log_info=self.log_info)
  #} end def

  def __del__(self): #{
    if (hasattr(self, "job_collections")): #{
      for job_collection in self.job_collections: #{
        job_collection.CloseFullFile()
        for job in job_collection: #{
          job.Close()
        #} end for
      #} end for
    #} end if
    CloseLogFile(self)
  #} end def

  def Run(self): #{
    LogMsg(self, "creating parallel hybrid read-support jobs...")
    start = time.time()
    LogMsg(self, "Processing groups...")
    process_start = time.time()
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    for group in parser: #{
      ExtremeDebugMsg(self, "Processing group %i..." % group.id)
      ExtremeDebugMsg(self, group.FullDataString().rstrip())
      group_start = time.time()
      self.ProcessGroup(group)
      ExtremeDebugMsg(self, "Time spent processing group: %s\n" %
        TimeSpent(group_start) + "--------------------------\n")
    #} end for
    LogMsg(self, "Time spent processing groups: %s" %
      TimeSpent(process_start))
    raise HybridSupportError("TEMPORARILY blocking job submission for testing")
    for chromosome in self.job_collections.Chromosomes(): #{
      sort_start = time.time()
      job_collection = self.job_collections[chromosome]
      job_collection.SortFileByPosition()
      LogMsg(self, "Time spent sorting chromosome %s job data: %s" %
        (chromosome, TimeSpent(sort_start)))
      if (None == self.options.cluster_head): #{
        self.RunJobsLocally(job_collection)
      else:
        job_collection.SplitJobs()
        jobs_path = self.CreateJobsFile(job_collection)
        self.SubmitJobs(jobs_path, chromosome)
      #} end if
    #} end for
    LogMsg(self, "Time spent creating parallel hybrid read-support jobs: "
      "%s" % TimeSpent(start))
  #} end def

  def ProcessGroup(self, group): #{
    for candidate in group.members: #{
      # get contig_coords
      #raise HybridSupportError("not implemented: get contig coords")
      contig_coords = candidate.ContigEventRegion()
      # split the contig coordinates if they are wider than a read and
      #   the candidate has an internal gap
      if (candidate.gap and candidate.GapIsInternal() and
          self.options.read_length < contig_coords.Span()): #{
        left_side = CoordPairCls(contig_coords.min, contig_coords.min)
        right_side = CoordPairCls(contig_coords.max, contig_coords.max)
        contig_coords = [left_side, right_side]
      #} end if
      if (candidate.gap): #{
        candidate.SplitGapBlocks()
      #} end if
      for i in [0,1]: #{
        is_ctg_start = (0 == i)
        candidate.breakpoints[i].contig_region = ContigRegionCls(
          candidate.contig_info.id, contig_coords,
          candidate.alignments[i].is_pos_strand, is_ctg_start)
        candidate.breakpoints[i].genome_blocks = candidate.blocks[i]
      #} end for
      #candidate.breakpointA.contig_region = ContigRegionCls(
      #  candidate.contig_info.id, contig_coords,
      #  candidate.align_info_A.strand, "start")
      #candidate.breakpointB.contig_region = ContigRegionCls(
      #  candidate.contig_info.id, contig_coords,
      #  candidate.align_info_B.strand, "end")
      #candidate.breakpointA.genome_blocks = candidate.blocks_A
      #candidate.breakpointB.genome_blocks = candidate.blocks_B
      sorted_breakpoints = SortBreakpoints(breakpoints[0], breakpoints[1])
      records = (None, None)
      for i in [0,1]: #{
        records[i] = HybridRecordCls(candidate.group_id,
          candidate.candidate_id, i, sorted_breakpoints[i].chr,
          sorted_breakpoints[i].contig_region)
        # get genome blocks for record
        #raise HybridSupportError("not implemented: get genome blocks for record %i" % i)
        self.GetGenomicBlocks(records[i], sorted_breakpoints[i].genome_blocks,
          candidate.ContigOverlap())
        chrom = records[i].chromosome
        job_collection = self.job_collections[chrom]
        job_collection.AddRecord(chrom, records[i].ToString())
      #} end for
      #record_A = HybridRecordCls(candidate.group_id, candidate.candidate_id,
      #  "A", sorted_breakpoints[0].chr, sorted_breakpoints[0].contig_region)
      #record_B = HybridRecordCls(candidate.group_id, candidate.candidate_id,
      #  "B", sorted_breakpoints[1].chr, sorted_breakpoints[1].contig_region)
      # get genome blocks for record A
      #raise HybridSupportError("not implemented: get genome blocks for record A")
      #self.GetGenomicBlocks(record_A, sorted_breakpoints[0].genome_blocks, ctg_overlap)
      # get genome blocks for record B
      #raise HybridSupportError("not implemented: get genome blocks for record B")
      #self.GetGenomicBlocks(record_B, sorted_breakpoints[1].genome_blocks, ctg_overlap)
      #self.job_collections[record_A.chromosome].AddRecord(record_A.chromosome,
      #  record_A.ToString())
      #self.job_collections[record_B.chromosome].AddRecord(record_B.chromosome,
      #  record_B.ToString())
    #} end for
  #} end def

  def GetGenomicBlocks(self, record, align_blocks, ctg_overlap): #{
    # calculate the search region length as 5/4 the expected fragment length
    #   plus any contig overlap in the alignments
    region_len = (1.25 * float(self.options.frag_len)) + ctg_overlap
    # get the search blocks by cutting or extending the alignment blocks
    #   so that their total length is correct
    genome_blocks = CutOrExtendBlocks(align_blocks, region_len,
      record.IsLeftBreakpoint())
    # merge any adjacent blocks
    genome_blocks = MergeAdjacentBlocks(genome_blocks)
    # add the search blocks to the record
    for genome_block in genome_blocks: #{
      record.AddGenomeBlock(genome_block)
    #} end for
  #} end def

  def RunJobsLocally(self, job_collection): #{
    LogMsg(self, "Running chromosome %s jobs locally..." %
      job_collection.chromosome)
    for job in job_collection: #{
      job.Close()
      # create job-file in case cluster resubmission occurs
      job_file = FileBoxCls(job.CommandPath(), "w",
        "cannot create job command file")
      job_file.WriteLine(job.Command(local=False))
      job_file.Close()
      out_file = {'path':job.StdOutPath(), 'mode':"w",
        'fail_msg':"cannot create stdout file for job"}
      err_file = {'path':job.StdErrPath(), 'mode':"w",
        'fail_msg':"cannot create stderr file for job"}
      job_command = job.Command(local=True)
      DebugMsg(self, job_command)
      return_code = RunCommandFromString(job_command,
        stdout=out_file, stderr=err_file, dpt=self.options.dpt)
      if (0 != return_code and not self.log_info['debug']): #{
        LogMsg(self, job_command)
      #} end if
      if (0 > return_code): #{
        raise HybridSupportError("Hybrid read-support job command was "
          "terminated by signal %i" % return_code)
      elif (0 < return_code):
        raise HybridSupportError("Error running hybrid read-support job "
          "command: %i" % return_code)
      #} end if
    #} end for
  #} end def

  def CreateJobsFile(self, job_collection): #{
    jobs_path = os.path.join(self.options.jobs_dir,
      "hybrid_%s_jobs" % job_collection.chromosome)
    LogMsg(self, "Creating jobs file: %s" % jobs_path)
    jobs_file = FileBoxCls(jobs_path, "w", "cannot create jobs file")
    for job in job_collection: #{
      job.Close()
      jobs_file.WriteLine(job.Command(local=False))
    #} end for
    jobs_file.Close()
    return jobs_file.path
  #} end def

  def SubmitJobs(self, jobs_path): #{
    LogMsg(self, "Submitting hybrid read-support jobs...")
    submit_options = GetOptions(submit)
    submit.SetSubmitOptions(submit_options, "%s-hyb" % self.options.lib,
      jobs_path, False, self.options)
    if (submit.CheckPaths(submit_options)): #{
      submitter = submit.JobSubmitterCls(submit_options,
        log_info=self.log_info)
      submitter.Submit()
    else:
      raise HybridSupportError("cannot submit jobs: "
        "errors in submit options")
    #} end if
  #} end def
#} end class

# every candidate creates two hybrid records
#   (for either side of the breakpoint)
class HybridRecordCls: #{
  def __init__(self, group_id, candidate_id, region_id, chromosome,
      contig_region): #{
    self.group_id     = group_id
    self.candidate_id = candidate_id
    self.region_id    = region_id
    self.chromosome   = chromosome
    # the contig region has contig ID, coords,
    #   contig-to-genome strand, and contig_part (start or end)
    self.contig_region = contig_region
    # use blocks to create a list of genomic regions that the contig aligns to
    #   each "region" is a pair of coordinates
    #   (CoordPairCls from common.coord_pair)
    self.genome_blocks = list()
    # genome_coords = (min(genome_coords), max(genome_coords))
    self.genome_coords = CoordPairCls()
  #} end def

  def AddGenomeBlock(self, new_block): #{
    self.genome_blocks.append(new_block)
    self.genome_coords.Union(new_block)
  #} end def

  def IsLeftBreakpoint(self): #{
    return IsLeftBreakpoint(self.contig_region.is_ctg_start,
      self.contig_region.is_pos_strand)
  #} end def
#} end class

class ContigRegionCls: #{
  def __init__(self, ctg_id, coords, is_pos_strand, is_ctg_start): #{
    self.ctg_id = ctg_id
    self.coords = coords
    self.is_pos_strand = is_pos_strand
    self.is_ctg_start   = is_ctg_start
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class HybridSupportError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Create and submit parallel jobs for calculating "
    "hybrid read-support")
  args = [ "LIB", "BARNACLE_FILE", "PAIR_TO_GENOME_FILE",
    "READ_TO_CONTIG_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--min-overlap",
                    type="int", metavar="N",
                    help="Require that reads overlap breakpoints by at "
                         "least N bp. [default: %default]")
  parser.add_option("--records-per-job",
                    type="int", metavar="N",
                    help="Process N records in each job. [default=%default]")
  cluster_group = OptionGroup(parser, "Cluster Options")
  cluster_group.add_option("-c", "--cluster-head",
                    dest="cluster_head",
                    help="Cluster head node to submit to.")
  cluster_group.add_option("--mem", "--memory",
                    help="The memory requirement of the jobs "
                         "[default:\"%default\"].")
  cluster_group.add_option("--hostname",
                    help="The hostname(s) to submit to "
                         "[default:\"%default\"].")
  cluster_group.add_option("--queue",
                    help="The queue(s) to submit to [default:\"%default\"].")
  cluster_group.add_option("-w", "--wall-time",
                    metavar="H:MM:SS",
                    help="The maximum time to spend on the job.")
  cluster_group.add_option("--email",
                    help="E-mail status updates on submitted jobs to the "
                         "given email address")
  parser.add_option_group(cluster_group)
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  misc_group.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  misc_group.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(min_overlap=5,
                      records_per_job=1000,
                      mem="3G",
                      hostname="q*",
                      queue="all.q",
                      dpt=False,
                      force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  (input_dir, input_file_name) = os.path.split(options.barnacle_path)
  options.jobs_dir = os.path.join(os.path.dirname(input_dir),
    "cluster_hyb")
  CheckFilePath(options.p2g_path, "pair-to-genome alignments", path_errors)
  CheckFilePath(options.r2c_path, "read-to-contigs alignments", path_errors)
  if (0 == len(path_errors)): #{
    # attempt to create output directory, if it does not already exist
    CheckDirPath(options.jobs_dir, "jobs directory", path_errors,
      create=True)
  #} end if
  options.log_file_name = GetLogPath(options.barnacle_path, "hyb.submit",
    options.jobs_dir)
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
    options.p2g_path      = EnsureAbsPath(args[2])
    options.r2c_path      = EnsureAbsPath(args[3])
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = HybridSubmitterCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.Run()
      except (MyError), e:
        ErrMsg("ERROR while creating parallel hybrid read-support jobs:\n"
          "  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "BARNACLE prediction data file for that library (BARNACLE_FILE); "
      "the path to a bam-formatted paired-read to genome alignment file "
      "(PAIR_TO_GENOME_FILE); and the path to a bam-formatted read "
      "to contig alignment file (READ_TO_CONTIG_FILE).")
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
