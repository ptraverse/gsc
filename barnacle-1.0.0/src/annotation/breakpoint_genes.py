#! /usr/bin/env python
"""
breakpoint_genes.py

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
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  ShouldChromUseChr, RunOverlapCode, NormalizeChrID)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls, CleanLine)
from utils.subprocesses import STREAM_OUT
from alignment_processing.gene_overlap import FeatureOverlapCls
from parsers.candidate_group_parser import CandidateGroupParserCls
from parsers.tokenizer import TokenizerCls
from common.grouped_candidate import ParseFullEventID

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "BREAKPOINT GENES SUCCESS"
MSG_FAIL = "BREAKPOINT GENES FAIL"

class BreakpointGeneAnnotatorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    self.overlaps_stream = None
    self.curr_overlap = None
    self.output_file  = None
  #} end def

  def __del__(self): #{
    if (hasattr(self, "output_file") and None != self.output_file): #{
      self.output_file.Close()
    #} end if
    CloseLogFile(self)
  #} end def

  def Run(self): #{
    LogMsg(self, "Adding breakpoint gene annotation...")
    start = time.time()
    group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    self.CreateGroupCoordsFile(group_parser)
    self.overlaps_stream = RunOverlapCode(self.options.genes_path,
      self.options.group_coords_path, STREAM_OUT, dpt=self.options.dpt,
      log_info=self.log_info)
    self.output_file = FileBoxCls(self.options.out_path, "w",
      "cannot create output file")
    overlaps_remain = True
    for group in group_parser: #{
      # get the breakpoint gene annotation for the group
      if (overlaps_remain): #{
        try: #{
          self.GetBreakpointOverlaps(group)
        except StopIteration:
          DebugMsg(self, "No overlaps remain")
          overlaps_remain = False
        #} end try
      #} end if
      # write the group to the output file
      self.output_file.Write(group.FullDataString())
    #} end for
    self.output_file.Close()
    if (not self.options.keep_coords_file): #{
      os.remove(self.options.group_coords_path)
    #} end def
    LogMsg(self, "Time spent adding breakpoint gene annotation: %s" %
      TimeSpent(start))
  #} end def

  def CreateGroupCoordsFile(self, group_parser): #{
    if (self.options.use_existing_group_coords): #{
      LogMsg(self, "Using existing group coordinates file.")
      return
    #} end if
    # check whether to use "chr" in chromosome names in coordinates file
    use_chr = ShouldChromUseChr(1, self.options.genes_path,
      "exon coordinates", self.log_info)
    # open the group coordinates file
    group_coords_file = FileBoxCls(self.options.group_coords_path, "w",
      "cannot create event coordinates file")
    for group in group_parser: #{
      ExtremeDebugMsg(self, "Writing coordinates for group %i" % group.id)
      self.WriteGroupCoords(group, group_coords_file, use_chr)
    #} end for
    group_parser.Close()
    group_coords_file.Close()
  #} end def

  def WriteGroupCoords(self, group, group_coords_file, use_chr): #{
    for candidate in group.members: #{
      if (candidate.gap): #{
        # write gap group coordinates
        self.WriteGapGroupCoords(candidate, group_coords_file, use_chr)
      else:
        # write split group coordinates
        self.WriteSplitGroupCoords(candidate, group_coords_file, use_chr)
      #} end if
    #} end for
  #} end def

  def WriteGapGroupCoords(self, candidate, group_coords_file, use_chr): #{
    gap_coords = ConstructBEDString(candidate.align_info_B.chrom, use_chr,
      candidate.align_info_B.genome_start, candidate.align_info_B.genome_end,
      "%sA" % candidate.IDString())
    group_coords_file.WriteLine(gap_coords)
  #} end def

  def WriteSplitGroupCoords(self, candidate, group_coords_file, use_chr): #{
    #split_coords_A = GroupCoordsCls(candidate.align_info_A.chrom,
    #  candidate.align_info_A.genome_end - self.options.event_buffer,
    #  candidate.align_info_A.genome_end + self.options.event_buffer,
    #  "%sA" % candidate.IDString(), use_chr)
    #group_coords_file.WriteLine("%s" % split_coords_A.ToString())
    #split_coords_B = GroupCoordsCls(candidate.align_info_B.chrom,
    #  candidate.align_info_B.genome_start - self.options.event_buffer,
    #  candidate.align_info_B.genome_start + self.options.event_buffer,
    #  "%sB" % candidate.IDString(), use_chr)
    #group_coords_file.WriteLine("%s" % split_coords_B.ToString())
    region_ids = ("A", "B")
    for i in [0,1]: #{
      genome_coords = (candidate.alignments[i].genome_start,
        candidate.alignments[i].genome_end)
      split_coords = ConstructBEDString(candidate.alignments[i].chrom, use_chr,
        genome_coords[1-i] - self.options.event_buffer,
        genome_coords[1-i] + self.options.event_buffer,
        "%s%s" % (candidate.IDString(), region_ids[i]))
      group_coords_file.WriteLine(split_coords)
    #} end for
  #} end def

  def GetBreakpointOverlaps(self, group): #{
    # clear any previous breakpoint genes
    group.ClearBPGenes()
    # skip overlaps for groups that come before the current group
    while (not hasattr(self, "curr_overlap") or None == self.curr_overlap or
        self.curr_overlap.group_id < group.id): #{
      self.GetNextOverlap()
    #} end while
    # create a dictionary of the members of the current group
    candidates_dict = dict((candidate.candidate_id, candidate) for
      candidate in group.members)
    # get all overlaps for the current group
    while (self.curr_overlap.group_id == group.id): #{
      ExtremeDebugMsg(self, "Found overlap for %i%s" %
        (self.curr_overlap.group_id, self.curr_overlap.candidate_id))
      if (self.curr_overlap.candidate_id in candidates_dict): #{
        ExtremeDebugMsg(self, "  candidate ID in dictionary!")
        self.AddBreakPointGene(candidates_dict[self.curr_overlap.candidate_id])
      #} end if
      self.GetNextOverlap()
    #} end while
  #} end def

  def GetNextOverlap(self): #{
    if (None == self.overlaps_stream): #{
      raise BreakpointGeneAnnotatorError("breakpoint gene overlap stream "
        "is not open!")
    #} end if
    overlap_line = CleanLine(self.overlaps_stream.next())
    tokenizer = TokenizerCls(overlap_line, delimiter=" ",
      log_info=self.log_info)
    try:
      self.curr_overlap = FeatureOverlapCls(tokenizer, multi_target=True)
    except ValueError,e:
      raise BreakpointGeneAnnotatorError("error parsing overlap line: "
        "%s\n%s" % (overlap_line, e))
    # end try
    ParseFullEventID(self.curr_overlap, self.curr_overlap.query_id)
  #} end def

  def AddBreakPointGene(self, candidate): #{
    ExtremeDebugMsg(self, "Adding breakpoint gene to candidate")
    if (None == self.curr_overlap): #{
      ExtremeDebugMsg(self, "  NO OVERLAP!")
      return
    #} end if
    if (self.curr_overlap.group_id != candidate.group_id): #{
      raise BreakpointGeneAnnotatorError("Group ID: %i does not match overlap "
        "ID: %i" % (candidate.group_id, self.curr_overlap.group_id))
    #} end if
    if (self.curr_overlap.candidate_id != candidate.candidate_id): #{
      raise BreakpointGeneAnnotatorError("Candidate ID: %s does not match "
        "overlap ID: %s" % (candidate.candidate_id,
        self.curr_overlap.candidate_id))
    #} end if
    ExtremeDebugMsg(self, "  OVERLAP: %s" % "\n  ".join((":".join((
      id_dict['gene'], id_dict['transcript'], id_dict['feature'],
      id_dict['num'])) for id_dict in self.curr_overlap.target_list)))
    candidate.AddGenes("breakpoint_%s" % self.curr_overlap.region_id,
      self.curr_overlap.target_list)
    ExtremeDebugMsg(self, "  AFTER: %s" % candidate.BPGenesString())
  #} end def
#} end class

def ConstructBEDString(chrom, use_chr, coord1, coord2, id): #{
  #if (use_chr): #{
  #  chrom = "chr%s" % chrom
  #else:
  #  chrom = chrom.replace("chr","")
  #} end if
  chrom = NormalizeChrID(chrom, use_chr)
  (left, right) = sorted((coord1, coord2))
  data_list = list([chrom, "%i" % left, "%i" % right, id, ])
  return " ".join(data_list)
#} end def

#### EXCEPTION CLASSES ####
class BreakpointGeneAnnotatorError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("") #TODO
  args = [ "LIB", "BARNACLE_FILE", "GENES_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--event-buffer",
                    type="int", metavar="N",
                    help="For split events, check whether the breakpoint "
                         "coordinates are within N bases of an exon. "
                         "[default: %default]")
  parser.add_option("--use-existing-coords",
                    action="store_true", dest="use_existing_group_coords",
                    help="If a breakpoint coordinates file already exists, "
                         "just use it rather than generating a new one.")
  parser.add_option("--keep-coords-file",
                    action="store_true",
                    help="After exon-overlap processing, do not remove the "
                         "alignment-coordinates file that was produced.")
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.set_defaults(event_buffer=5,
                      use_existing_group_coords=False,
                      keep_coords_file=False,
                      dpt=False,
                      force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (0 > options.event_buffer): #{
    ErrMsg("Invalid event buffer value: \"%i\": must be non-negative." %
      options.event_buffer)
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle input data", path_errors)
  CheckFilePath(options.genes_path, "gene feature annotation", path_errors)
  # get the input directory
  (input_dir, input_file_name) = os.path.split(options.barnacle_path)
  # get the output path
  output_dir = GetOutDir(input_dir, "breakpoint_genes")
  options.out_path = os.path.join(output_dir, input_file_name)
  # get the path for the group coordinates and overlaps output files
  input_root = os.path.splitext(input_file_name)[0]
  options.group_coords_path = os.path.join(output_dir,
    "%s.coords.bed" % input_root)
  if (options.use_existing_group_coords and
      not os.path.isfile(options.group_coords_path)): #{
    options.use_existing_group_coords = False
  #} end if
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "breakpoint_genes", output_dir)
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
    options.genes_path    = EnsureAbsPath(args[2])
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = BreakpointGeneAnnotatorCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.Run()
      except (MyError), e:
        ErrMsg("ERROR while adding breakpoint gene annotation:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
     "Barnacle data file for that library (BARNACLE_FILE); and the path "
     "to a BED-format gene feature coordinates file (GENES_FILE).")
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
