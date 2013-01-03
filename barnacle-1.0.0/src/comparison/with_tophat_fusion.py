#! /usr/bin/env python
"""
with_tophat_fusion.py

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
from utils.general import (SetupMainClass, TimeSpent, WriteCommand)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls
from common.breakpoint import BreakpointCls, SortBreakpoints

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

class TopHatComparisonCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    # barnacle_dict[chrA][chrB][coord_keyA][coord_keyB] =
    #   list of Barnacle events
    # REMINDER: USE MULTI-DICT FUNCTIONS
    self.barnacle_dict = dict()
    # a simple list, for efficiently outputting Barnacle-only events
    self.barnacle_list = list()
    self.counts = {'barnacle':0, 'tophat':0, 'common':0}
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def CompareResults(self): #{
    LogMsg(self, "Comparing Barnacle and TopHat-Fusion results...")
    start = time.time()
    # load the Barnacle results
    self.LoadBarnacleResults()
    # open the output files (barnacle_only, tophatfusion_only, common)
    self.OpenOutputFiles()
    for tophat_path in self.options.tophat_paths: #{
      # open TopHat-Fusion results file
      tophat_file = TopHatFileCls(tophat_path, log_info=self.log_info)
      for tophat_event in tophat_file: #{
        DebugMsg(self, "Processing TopHat-Fusion event...\n  %s" %
          tophat_event.ToString())
        self.ProcessTopHatEvent(tophat_event)
      #} end for
    #} end for
    # write Barnacle-only events
    self.WriteBarnacleOnlyEvents()
    LogMsg(self, "Time spent comparing: %s" % TimeSpent(start))
    LogMsg(self, "Found %i Barnacle specific events" %
      self.counts['barnacle'])
    LogMsg(self, "Found %i TopHat-Fusion specific events" %
      self.counts['tophat'])
    LogMsg(self, "Found %i common events" % self.counts['common'])
  #} end def

  def LoadBarnacleResults(self): #{
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    for group in parser: #{
      # mark the group as Barnacle-only by default
      group.barnacle_only = True
      self.SetupBarnacleEventBreakpoints(group)
      chr1 = group.breakpoint1.chr
      # update the dictionary
      if (chr1 not in self.barnacle_dict): #{
        self.barnacle_dict[chr1] = dict()
      #} end if
      chr2 = group.breakpoint2.chr
      if (chr2 not in self.barnacle_dict[chr1]): #{
        self.barnacle_dict[chr1][chr2] = dict()
      #} end if
      coord_key1 = group.breakpoint1.coord / self.options.coord_buffer
      if (coord_key1 not in self.barnacle_dict[chr1][chr2]): #{
        self.barnacle_dict[chr1][chr2][coord_key1] = dict()
      #} end if
      coord_key2 = group.breakpoint2.coord / self.options.coord_buffer
      if (coord_key2 not in self.barnacle_dict[chr1][chr2][coord_key1]): #{
        self.barnacle_dict[chr1][chr2][coord_key1][coord_key2] = list()
      #} end if
      DebugMsg(self, "Adding Barnacle event (chr1:%s, chr2:%s, " %
        (chr1, chr2) + "key1:%i, key2:%i)" % (coord_key1, coord_key2))
      self.barnacle_dict[chr1][chr2][coord_key1][coord_key2].append(group)
      self.barnacle_list.append(group)
    #} end for
  #} end def

  def SetupBarnacleEventBreakpoints(self, group): #{
    # use the first member
    member = group.members[0]
    if (hasattr(member, 'breakpoint_genes_A')): #{
      gene_set = member.breakpoint_genes_A
    else:
      gene_set = member.genes_A
    #} end if
    member.breakpointA.gene = member.GeneSetString(gene_set)
    member.breakpointA.end = member.ends[0]
    if (hasattr(member, 'breakpoint_genes_B')): #{
      gene_set = member.breakpoint_genes_B
    else:
      gene_set = member.genes_B
    #} end if
    member.breakpointB.gene = member.GeneSetString(gene_set)
    member.breakpointB.end = member.ends[1]
    (group.breakpoint1, group.breakpoint2) = (
      SortBreakpoints(member.breakpointA, member.breakpointB))
  #} end def

  def OpenOutputFiles(self): #{
    path = os.path.join(self.options.output_dir, "barnacle_only.txt")
    self.barnacle_out_file = FileBoxCls(path, "w",
      "cannot create Barnacle only results file")
    path = os.path.join(self.options.output_dir, "tophat_only.txt")
    self.tophat_out_file = FileBoxCls(path, "w",
      "cannot create TopHat-Fusion only results file")
    path = os.path.join(self.options.output_dir, "common.txt")
    self.common_out_file = FileBoxCls(path, "w",
      "cannot create common results file")
  #} end def

  def ProcessTopHatEvent(self, tophat_event): #{
    if (self.TopHatEventInBarnacle(tophat_event)): #{
      DebugMsg(self, "Common event detected!")
      self.common_out_file.WriteLine("TopHat-Fusion: %s" %
        tophat_event.ToString())
      self.common_out_file.WriteLine("-"*80)
      self.counts['common'] += 1
    else:
      DebugMsg(self, "TopHat-Fusion specific event!")
      self.tophat_out_file.WriteLine(tophat_event.ToString())
      self.counts['tophat'] += 1
    #} end if
  #} end def

  def TopHatEventInBarnacle(self, tophat_event): #{
    if (tophat_event.breakpoint1.chr not in self.barnacle_dict): #{
      DebugMsg(self, "TopHat event chr1 (%s) not in Barnacle events" %
        tophat_event.breakpoint1.chr)
      return False
    #} end if
    chr1_dict = self.barnacle_dict[tophat_event.breakpoint1.chr]
    if (tophat_event.breakpoint2.chr not in chr1_dict): #{
      DebugMsg(self, "TopHat event chr2 (%s) not in Barnacle events" %
        tophat_event.breakpoint2.chr)
      return False
    #} end if
    chr2_dict = chr1_dict[tophat_event.breakpoint2.chr]
    event_in_barnacle = False
    coord_key1 = tophat_event.breakpoint1.coord / self.options.coord_buffer
    # also check the bins before and after
    for key1 in range(coord_key1-1,coord_key1+2): #{
      if (key1 not in chr2_dict): #{
        DebugMsg(self, "TopHat event coord1-key (%i) not in Barnacle events" %
          key1)
        continue
      #} end if
      coord1_dict = chr2_dict[key1]
      coord_key2 = tophat_event.breakpoint2.coord / self.options.coord_buffer
      # also check the bins before and after
      for key2 in range(coord_key2-1,coord_key2+2): #{
        if (key2 not in coord1_dict): #{
          DebugMsg(self, "TopHat event coord2-key (%i) " % key2 +
            "not in Barnacle events")
          continue
        #} end if
        coord2_list = coord1_dict[key2]
        for barnacle_event in coord2_list: #{
          if (self.CompareEvents(barnacle_event, tophat_event)): #{
            event_in_barnacle = True
            barnacle_event.barnacle_only = False
            self.WriteCommonBarnacleEvent(barnacle_event)
          #} end if
        #} end for
      #} end for
    #} end for
    return event_in_barnacle
  #} end def

  def CompareEvents(self, barnacle_event, tophat_event): #{
    DebugMsg(self, "Comparing TopHat-Fusion event to Barnacle event %i" %
      barnacle_event.id)
    dist1 = abs(barnacle_event.breakpoint1.coord - tophat_event.breakpoint1.coord)
    DebugMsg(self, "Distance for breakpoint1: %i" % dist1)
    if (dist1 > self.options.coord_buffer): #{
      return false
    #} end if
    dist2 = abs(barnacle_event.breakpoint2.coord - tophat_event.breakpoint2.coord)
    DebugMsg(self, "Distance for breakpoint2: %i" % dist2)
    if (dist2 > self.options.coord_buffer): #{
      return false
    #} end if
    return True
  #} end def

  def WriteCommonBarnacleEvent(self, barnacle_event): #{
    gene_strs = [
      "%s(%s')" % (barnacle_event.breakpoint1.gene,
        barnacle_event.breakpoint1.end),
      "%s(%s')" % (barnacle_event.breakpoint2.gene,
        barnacle_event.breakpoint2.end),
    ]
    coord_strs = [
      "%s:%i" % (barnacle_event.breakpoint1.chr,
        barnacle_event.breakpoint1.coord),
      "%s:%i" % (barnacle_event.breakpoint2.chr,
        barnacle_event.breakpoint2.coord),
    ]
    barnacle_str = ("G%i %s %s" % (barnacle_event.id,
      ";".join(gene_strs), ";".join(coord_strs)))
    self.common_out_file.WriteLine("Barnacle: %s" % barnacle_str)
  #} end def

  def WriteBarnacleOnlyEvents(self): #{
    for group in self.barnacle_list: #{
      if (group.barnacle_only): #{
        self.barnacle_out_file.Write(group.FullDataString())
        self.counts['barnacle'] += 1
      #} end if
    #} end for
  #} end def
#} end class

class TopHatFileCls: #{
  def __init__(self, path, log_info=None): #{
    self.file = FileBoxCls(path, "r", "cannot read TopHat-Fusion results file")
    self.log_info = log_info
  #} end def

  def __del__(self): #{
    self.close()
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  def next(self): #{
    # the first line should start with "allAtOnce" and
    # it contains the breakpoint coordinates
    # parse the tophat line
    tophat_event = TopHatEventCls(self.file.next())
    # the next two lines should be "sequence" lines
    tophat_event.CheckSeqLine(self.file.next())
    tophat_event.CheckSeqLine(self.file.next())
    # the next lines should be... scores?
    tophat_event.CheckScoreLine(self.file.next())
    # the next line should have the gene ids
    tophat_event.ParseGenesLine(self.file.next())
    # skip the final line
    self.file.next()
    return tophat_event
  #} end def

  def close(self): #{
    if (hasattr(self, "file") and None != self.file and
        not self.file.closed): #{
      self.file.close()
    #} end if
  #} end def
#} end class

class TopHatEventCls: #{
  def __init__(self, start_line): #{
    if (not start_line.startswith("allAtOnce")): #{
      raise TopHatComparisonError("First line of TopHat-Fusion result should "
        "start with \"allAtOnce\": %s" % start_line)
    #} end if
    fields = start_line.split(" ")
    chroms = fields[1]
    (chrA, chrB) = chroms.split("-")
    coordA = int(fields[2])
    breakpointA_str = "%s:%i(up)" % (chrA, coordA)
    self.breakpointA = BreakpointCls(breakpointA_str)
    coordB = int(fields[3])
    breakpointB_str = "%s:%i(up)" % (chrB, coordB)
    self.breakpointB = BreakpointCls(breakpointB_str)
    self.seq_patt = re.compile(r"^[ACGT]{50,50} [ACGT]{50,50}$")
    self.score_patt = re.compile(r"^[0-9]{50,50} [0-9]{50,50}$")
  #} end def

  def CheckSeqLine(self, seq_line): #{
    if (None == self.seq_patt.search(seq_line)): #{
      raise TopHatComparisonError("Invalid sequence line: %s" % seq_line)
    #} end if
  #} end def

  def CheckScoreLine(self, score_line): #{
    if (None == self.score_patt.search(score_line)): #{
      raise TopHatComparisonError("Invalid score line: %s" % score_line)
    #} end if
  #} end def

  def ParseGenesLine(self, genes_line): #{
    fields = genes_line.split(" ")
    self.breakpointA.gene = fields[0]
    self.breakpointA.part = fields[1].split("(")[0]
    self.breakpointB.gene = fields[2]
    self.breakpointB.part = fields[3].split("(")[0]
    (self.breakpoint1, self.breakpoint2) = (
      SortBreakpoints(self.breakpointA, self.breakpointB))
  #} end def

  def ToString(self): #{
    gene_strs = list()
    coord_strs = list()
    if (hasattr(self, "breakpoint1") and hasattr(self, "breakpoint2")): #{
      gene_strs.extend([
        "%s(%s)" % (self.breakpoint1.gene, self.breakpoint1.part),
        "%s(%s)" % (self.breakpoint2.gene, self.breakpoint2.part),
      ])
      coord_strs.extend([
        "%s:%i" % (self.breakpoint1.chr, self.breakpoint1.coord),
        "%s:%i" % (self.breakpoint2.chr, self.breakpoint2.coord),
      ])
    else:
      gene_strs.extend(["N/A", "N/A"])
      coord_strs.extend([
        "%s:%i" % (self.breakpointA.chr, self.breakpointA.coord),
        "%s:%i" % (self.breakpointB.chr, self.breakpointB.coord),
      ])
    #} end if
    return "%s %s" % ("/".join(gene_strs), "/".join(coord_strs))
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class TopHatComparisonError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Compares a list of fusions predicted by Barnacle "
    "and a list predicted by TopHat-Fusion and returns which predictions "
    "are common to both, or unique to one or the other.")
  args = [ "LIB", "BARNACLE_FILE", "TOPHAT_FUSION_FILES", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--coord-buffer",
                    type="int", metavar="N",
                    help="Consider events identical if their coordinates "
                         "are within Nbp of each other. [default: %default]")
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(coord_buffer=1000,
                      force=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (not opts_good): #{
    ErrMsg("bad option") #TODO
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  for tophat_path in options.tophat_paths: #{
    CheckFilePath(tophat_path, "TopHat-Fusion results", path_errors)
  #} end for
  # get and check the output path
  input_dir = os.path.dirname(options.barnacle_path)
  options.output_dir = os.path.join(input_dir, "tophat_comparison")
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "tophat_comparison", options.output_dir) #TODO
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
    options.tophat_paths = [EnsureAbsPath(path) for
      path in args[2].strip(",").split(",")]
    if (CheckPaths(options)): #{
      try: #{
        event_comparer = TopHatComparisonCls(options)
        WriteCommand(event_comparer, sys.argv)
        event_comparer.CompareResults()
      except (MyError), e:
        ErrMsg("ERROR while comparing results:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle data file for that library (BARNACLE_FILE); and a "
      "comma-delimited list of paths to TopHat-Fusion results files "
      "(TOPHAT_FUSION_FILES).")
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
