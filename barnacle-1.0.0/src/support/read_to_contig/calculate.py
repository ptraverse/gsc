#! /usr/bin/env python
"""
calculate.py

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
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand, GetGroupID,
  CheckConfigCommands)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  FileBoxCls, CleanLine)
from utils.subprocesses import SortFile
from parsers.tokenizer import GetFieldValue
from support.samtools import SAMToolsCls
from support.SAM_record import SAM_record
from group import R2CGroupCls
from align import R2CAlignCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "R2C CALCULATE SUCCESS"
MSG_FAIL = "R2C CALCULATE FAIL"

class R2CSupportCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    CheckConfigCommands(self, "samtools")
    self.groups = dict()
    self.contigs = dict()
    self.min_missing_alt = None
    self.num_dup_aligns = 0
    self.num_multi_aligns = 0
    self.found_reads = False
    self.num_inconsistent1 = 0
    self.num_inconsistent2 = 0
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def Run(self): #{
    LogMsg(self, "Calculating read-to-contig support...")
    start = time.time()
    # get groups, contigs, and regions from data file
    self.ParseGroupsFile()
    # for each contig
    LogMsg(self, "Processing contigs...")
    process_start = time.time()
    for contig in self.contigs.itervalues(): #{
      # find all the reads supporting that contig
      self.ProcessContig(contig)
    #} end for
    LogMsg(self, "Time spent processing contigs: %s" %
      TimeSpent(process_start))
    # create output file
    self.CreateOuputFile()
    # for each group
    num_members_with_support = 0
    LogMsg(self, "Writing output...")
    write_start = time.time()
    for group_id in sorted(self.groups.keys()): #{
      # output the support for each member
      # (number of reads supporting each window)
      for member in self.groups[group_id].members.itervalues(): #{
        # Cannot calculate min and average support until integrating
        #member.CalculateSupport()
        if (member.HasSupport()): #{
          num_members_with_support += 1
          self.output_file.WriteLine(member.SupportString())
        #} end if
      #} end for
    #} end for
    LogMsg(self, "Time spent writing output: %s" %
      TimeSpent(write_start))
    # sort the results by the member id
    #self.SortOutputByMember()
    LogMsg(self, "Time spent calculating read-to-contig support: %s" %
      TimeSpent(start))
    LogMsg(self, "Found support for %i members" % num_members_with_support)
    if (None != self.min_missing_alt): #{
      LogMsg(self, "Warning: alternate hits missing for reads with %i " %
        self.min_missing_alt + "or more hits.")
    #} end if
    if (0 < self.num_dup_aligns): #{
      LogMsg(self, "Warning: %i reads have the primary alignment duplicated in "
        "their alternative hits." % self.num_dup_aligns)
    #} end if
    if (0 < self.num_multi_aligns): #{
      LogMsg(self, "Warning: %i reads have multiple alignments to the same "
        "contig (only the first encountered will be used)." %
        self.num_multi_aligns)
    #} end if
    if (0 < (self.num_inconsistent1 + self.num_inconsistent2)): #{
      LogMsg(self, "Warning: %i (%i) reads have inconsistent hit counting." %
          (self.num_inconsistent1, self.num_inconsistent2))
    #} end if
    if (0 == num_members_with_support): #{
      LogMsg(self, "WARNING: Did not find support for any members!")
    #} end if
    if (not self.found_reads): #{
      raise R2CSupportError("Could not find reads aligned to any of the "
        "contigs being considered!")
    #} end if
    # remove the temporary samtools output files
    temp_sam_path = os.path.join(self.options.output_dir, "sam_out_tmp.err")
    if (os.path.isfile(temp_sam_path)): #{
      os.remove(temp_sam_path)
    #} end if
  #} end def

  def ParseGroupsFile(self): #{
    LogMsg(self, "Parsing groups file...")
    start = time.time()
    fail_msg = "cannot open groups file"
    groups_file = FileBoxCls(self.options.barnacle_path, "r", fail_msg)
    for member_line in groups_file: #{
      ExtremeDebugMsg(self, "%s\nMEMBER LINE: %s" % ("-"*20, member_line))
      (ctg_id, member_id, region_coords, mates_str) = member_line.split(" ")
      group_id = GetGroupID(member_id)
      mates = GetFieldValue(mates_str)
      if (group_id not in self.groups): #{
        new_group = R2CGroupCls(group_id, log_info=self.log_info)
        new_group.AddMember(member_id, ctg_id, region_coords, mates)
        new_group.InitializeWindows(self.options.min_overlap)
        ExtremeDebugMsg(self, "GROUP ADDED:\n%s" % new_group.DebugString())
        self.groups[group_id] = new_group
      #} end if
      self.UpdateContigs(ctg_id, member_id, region_coords, primary=True)
      self.UpdateSecondaryContigs(mates)
    #} end for
    ExtremeDebugMsg(self, "-"*20)
    LogMsg(self, "Time spent parsing groups file: %s" % TimeSpent(start))
  #} end def

  def UpdateContigs(self, ctg_id, member_id, region_coords, primary): #{
    if (ctg_id not in self.contigs): #{
      self.contigs[ctg_id] = R2CContigCls(ctg_id, member_id, region_coords,
        primary, self.options, log_info=self.log_info)
      ExtremeDebugMsg(self, "CONTIG ADDED: %s" %
        self.contigs[ctg_id].DebugString())
    else:
      self.contigs[ctg_id].Update(member_id, region_coords, primary=primary)
      ExtremeDebugMsg(self, "CONTIG UPDATED: %s" %
        self.contigs[ctg_id].DebugString())
    #} end if
  #} end def

  def UpdateSecondaryContigs(self, mates): #{
    if ("none" == mates.lower()): #{
      return
    #} end if
    for mate_str in mates.split(";"): #{
      (ctg_id, member_id, region_coords) = mate_str.split(",")
      self.UpdateContigs(ctg_id, member_id, region_coords, primary=False)
    #} end for
  #} end def

  def ProcessContig(self, contig): #{
    ExtremeDebugMsg(self, "Processing contig %s" % contig.ctg_id)
    if (not contig.primary): #{
      ExtremeDebugMsg(self, "  Skipping secondary contig.")
      return
    #} end if
    # for each read mapping to the current contig
    #for read in contig.GetReads(self.options.r2c_path, self.options): #{
    for read in contig.GetReads(): #{
      self.found_reads = True
      # check which members read supports
      self.ProcessRead(read, contig)
      self.ProcessAlternativeHits(read)
    #} end for
  #} end def

  def ProcessRead(self, read, contig): #{
    ExtremeDebugMsg(self, "  Processing read %s" % read.id)
    unique_by_group = dict()
    for group_id in contig.groups: #{
      if (group_id in self.groups): #{
        ctgs_in_group = self.groups[group_id].contigs
        unique_by_group[group_id] = read.IsUnique(ctgs_in_group)
      #} end if
    #} end for
    for member_id in contig.members: #{
      self.ProcessMember(read, member_id, unique_by_group)
    #} end for
  #} end def

  def ProcessMember(self, read, member_id, unique_by_group): #{
    group_id = GetGroupID(member_id)
    member = self.groups[group_id].members[member_id]
    ExtremeDebugMsg(self, "    Processing member %s: %i-%i" % (member_id,
      member.left, member.right))
    unique = unique_by_group[group_id]
    member.ProcessRead(read, unique)
  #} end def

  def ProcessAlternativeHits(self, read): #{
    ExtremeDebugMsg(self, "  Processing alternative hits for read %s" %
      read.id)
    if (1 == read.num_best_hits or read.is_alt_hit): #{
      ExtremeDebugMsg(self, "  No alternative hits")
      return
    elif (read.alt_missing):
      ExtremeDebugMsg(self, "    Alternative hits missing for read %s "
        "with %i hits" % (read.id, read.num_best_hits))
      if (None == self.min_missing_alt or
          read.num_best_hits < self.min_missing_alt): #{
        self.min_missing_alt = read.num_best_hits
      #} end if
      return
    #} end if
    if ((read.num_best_hits+read.num_subopt) != len(read.alt_hits)+1): #{
      #raise R2CSupportError("inconsistent hit counting for read %s:\n  " %
      self.num_inconsistent1 += 1
      ExtremeDebugMsg(self, "Warning: inconsistent hit counting for read %s:" %
        read.id + "\n  num best hits: %i " % read.num_best_hits +
        "num subopt hits: %i " % read.num_subopt +
        "alternate hits: %s" % ";".join(read.alt_hits))
    #} end if
    num_alt_best_hits = 0
    read.curr_alt = 0
    while (read.curr_alt < len(read.alt_hits)): #{
      new_read = R2CAlignCls(read)
      read.curr_alt += 1
      if (new_read.edit_dist > read.edit_dist):
        ExtremeDebugMsg(self, "    Suboptimal hit (%i > %i), skipping" %
          (new_read.edit_dist, read.edit_dist))
        continue
      # end if
      self.ProcessAlternativeHit(new_read, read)
      num_alt_best_hits += 1
    #} end while
    if (read.num_best_hits != num_alt_best_hits+1):
      #raise R2CSupportError("inconsistent hit counting for read %s:\n" %
      self.num_inconsistent2 += 1
      ExtremeDebugMsg(self, "Warning: inconsistent hit counting for read %s:" %
        read.id + "\n  num best hits: %i alternate best hits: %i" %
        (read.num_best_hits, num_alt_best_hits))
    # end if
  #} end def

  def ProcessAlternativeHit(self, new_read, old_read): #{
    if (ReadsAreEquivalent(new_read, old_read)): #{
      DebugMsg(self, "Warning: alignment for read %s to %s duplicated "
        "in alternative hit" % (new_read.id, new_read.ctg_id))
      ExtremeDebugMsg(self, "  %s" % old_read.ToString())
      self.num_dup_aligns += 1
      return
    elif (new_read.ctg_id == old_read.ctg_id):
      self.num_multi_aligns += 1
      ExtremeDebugMsg(self, "Warning: read %s " % new_read.id +
        "has multiple alignments to %s, " % new_read.ctg_id +
        "only the first encountered will be used.")
      return
      #raise R2CSupportError("Read %s already present in contig %s" %
      #  (new_read.id, new_read.ctg_id))
    #} end if
    if (not new_read.is_perfect): #{
      ExtremeDebugMsg(self, "    Skipping imperfect alignment")
      return
    #} end if
    if (new_read.ctg_id not in self.contigs): #{
      ExtremeDebugMsg(self, "    Contig %s not in groups, skipping" %
        new_read.ctg_id)
      return
    #} end if
    # process the alternative read
    contig = self.contigs[new_read.ctg_id]
    self.ProcessRead(new_read, contig)
  #} end def

  def CreateOuputFile(self): #{
    fail_msg = "cannot create read-to-contig support output file"
    self.output_file = FileBoxCls(self.options.output_path, "w", fail_msg)
  #} end def

  def SortOutputByMember(self): #{
    # ensure the output file is closed
    self.output_file.Close()
    sort_start = time.time()
    SortFile(self, self.output_file.path, 2, numeric=True)
    DebugMsg(self, "Time spent sorting output: %s" % TimeSpent(sort_start))
  #} end def
#} end class

class R2CContigCls: #{
  def __init__(self, ctg_id, member_id, region_coords, primary, options,
      log_info=None): #{
    self.log_info = log_info
    self.ctg_id  = ctg_id
    self.groups = set()
    self.groups.add(GetGroupID(member_id))
    self.members = set()
    self.members.add(member_id)
    (self.left, self.right) = map(int, region_coords.split("-"))
    # bool primary: is this a primary contig for the current job
    self.primary = primary
    self.samtools = SAMToolsCls(options.r2c_path, options,
      self.log_info, paired=False)
  #} end def

  def Update(self, member_id, region_coords, primary): #{
    self.groups.add(GetGroupID(member_id))
    self.members.add(member_id)
    (new_left, new_right) = map(int, region_coords.split("-"))
    self.left  = min(self.left,  new_left)
    self.right = max(self.right, new_right)
    if (primary): #{
      self.primary = True
    #} end if
  #} end def

  #def GetReads(self, r2c_path, ctgs_in_group):
  def GetReads(self): #{
    region = "%s:%i-%i" % (self.ctg_id, self.left, self.right)
    ExtremeDebugMsg(self, "Getting reads aligning to %s..." % region)
    for raw_record in self.samtools.GetReads(region): #{
      raw_record = CleanLine(raw_record)
      ExtremeDebugMsg(self, "%s" % raw_record)
      # use Rod's SAM_record code to parse reads from r2c alignment file
      record = SAM_record(raw_record)
      # make sure that the line was not a header (just in case)
      if (record.header): #{
        ExtremeDebugMsg(self, "Line \"%s\" is a header line, it will be "
          "skipped..." % record.raw_sam)
        continue
      #} end if
      # make sure the read is mapped (just in case)
      if (record.q_unmapped()): #{
        ExtremeDebugMsg(self, "Read %s is unmapped, will be skipped..." %
          record.qname)
        continue
      #} end if
      # make sure the read is a perfect match but allow soft-clipping
      if (not record.perfect_softclipped()): #{
        ExtremeDebugMsg(self, "Skipping imperfect alignment for read %s: %s "
          "ED:%s" % (record.qname, record.cigar, record.edit_distance))
        continue
      #} end if
      # use ctgs_in_group to mark reads as "unique" or not
      #new_read.InitializeFromSAMRecord(record, ctgs_in_group)
      # delay determining "unique" until group.CalculateReadToContigSupport()
      new_read = R2CAlignCls(record)
      ExtremeDebugMsg(self, new_read.ToString())
      yield new_read
    #} end for
  #} end def

  def DebugString(self): #{
    data_list = [
      self.ctg_id,
      ",".join(sorted(self.members)),
      "%i-%i" % (self.left, self.right),
      "P:%s" % self.primary,
    ]
    return " ".join(data_list)
  #} end def
#} end class

def ReadsAreEquivalent(read1, read2):
  return (read1.ctg_id    == read2.ctg_id and
          read1.strand    == read2.strand and
          read1.left      == read2.left   and
          read1.cigar     == read2.cigar  and
          read1.edit_dist == read2.edit_dist)
# end def

#### EXCEPTION CLASSES ####
class R2CSupportError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Calculates read-to-contig support for groups in "
    "the input file")
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
  parser.set_defaults(min_overlap=5,
                      dpt=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (not opts_good): #{
    ErrMsg("bad option") # TODO
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  CheckFilePath(options.r2c_path, "read-to-contigs alignments", path_errors)
  options.output_path = options.barnacle_path.replace(".data", ".out")
  options.output_dir = os.path.dirname(options.barnacle_path)
  CheckDirPath(options.output_dir, "output", path_errors, create=False)
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
        r2c_calculator = R2CSupportCls(options)
        WriteCommand(r2c_calculator, sys.argv)
        r2c_calculator.Run()
      except (MyError), e:
        ErrMsg("ERROR while calculating read-to-contig support:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle read-to-contig group data file (BARNACLE_FILE); and the "
      "path to a bam-formatted read to contig alignment file "
      "(READ_TO_CONTIG_FILE).")
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
