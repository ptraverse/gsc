#! /usr/bin/env python
"""
exon_bounds.py

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
  NormalizeChrID, NonStandardChr, OtherSide, SIDES)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.multi_dict import AddToMultiDict
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls
from parsers.genes.annotation import GeneAnnotationParserCls
from common.coord_pair import CoordPairCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "EXON BOUNDARY COUNTING SUCCESS"
MSG_FAIL = "EXON BOUNDARY COUNTING FAIL"

class ExonBoundCounterCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    # exon_bound_coords[chrom][prime_side][coord1][coord2] = gene_list
    self.exon_bound_coords = dict()
    self.output_file = None
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
    if (None != self.output_file): #{
      self.output_file.Close()
    #} end if
  #} end def

  def CountExonBoundaryBreakpoints(self): #{
    LogMsg(self, "Counting exon-boundary matching breakpoint coordinates...")
    start = time.time()
    # load the exon boundary coordinates from the annotations file(s)
    self.LoadExonBoundaryCoordinates()
    # create the parser
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    # create the output file
    self.CreateOutputFile()
    # iterate through the events
    LogMsg(self, "Processing events...")
    process_start = time.time()
    for event in parser: #{
      # count member breakpoint coordinates matching exon boundaries
      self.ProcessEvent(event)
    #} end for
    LogMsg(self, "Time spent processing events: %s\nTotal time spent "
      "counting exon-boundary matching breakpoint coordinates: %s" %
      (TimeSpent(process_start), TimeSpent(start)))
  #} end def

  def LoadExonBoundaryCoordinates(self): #{
    LogMsg(self, "Loading exon boundary coordinates...")
    start = time.time()
    # load the exon boundary coordinates from the main annotations file
    self.LoadExonBoundaryCoordinatesFromFile(self.options.genes_path)
    # load exon boundary coordinates from additional annotations files?
    LogMsg(self, "Time spent loading exon boundary coordinates: %s" %
      TimeSpent(start))
  #} end def

  def LoadExonBoundaryCoordinatesFromFile(self, path): #{
    DebugMsg(self, "Gene annotation file: %s" % path)
    # open the annotations file
    annotations_file = GeneAnnotationParserCls(path, log_info=self.log_info)
    skipped_chroms = set()
    # get the coordinates from the file
    for transcript in annotations_file: #{
      # fix chromosome names, if needed
      chrom = NormalizeChrID(transcript.chrom)
      if (NonStandardChr(chrom)): #{
        ExtremeDebugMsg(self, "Skipping transcript in strange chromosome: "
          "%s (%s)" % (chrom, transcript.chrom))
        skipped_chroms.add(chrom)
        continue
      #} end if
      prev_exon = None
      for (index, exon) in enumerate(transcript.SortedExons()): #{
        (exon.left, exon.right) = (exon.min, exon.max)
        # assume that exon list is sorted by left coordinate
        if (None != prev_exon and prev_exon.min > exon.min): #{
          raise ExonBoundCounterError("Transcript %s exons are not in "
            "order: %s, %s" % (transcript.transcript_id,
            prev_exon.ToString(), exon.ToString()))
        #} end if
        prev_exon = exon
        # do not include the left side of the first exon
        if (0 == index): #{
          exon.left = None
        #} end if
        # do not include the right side of the last exon
        if (len(transcript.exons) == (index+1)): #{
          exon.right = None
        #} end if
        # exon_bound_coords[chrom][prime_side][coord1][coord2] = gene_list
        for side in SIDES: #{
          if (None != getattr(exon, side)): #{
            keys = [chrom, side, getattr(exon, side),
              getattr(exon, OtherSide(side))]
            AddToMultiDict(self.exon_bound_coords, keys,
              transcript.transcript_id)
          #} end if
        #} end for
      #} end for
    #} end for
    if (0 < len(skipped_chroms)): #{
      DebugMsg(self, "Skipped transcripts in chromosomes: %s" %
        ", ".join(sorted(skipped_chroms)))
    #} end if
    # close the file
    annotations_file.close()
  #} end def

  def CreateOutputFile(self): #{
    input_file_name = os.path.basename(self.options.barnacle_path)
    output_path = os.path.join(self.options.output_dir, input_file_name)
    fail_msg = "could not create exon-boundary output file"
    self.output_file = FileBoxCls(output_path, "w", fail_msg)
  #} end def

  def ProcessEvent(self, event): #{
    DebugMsg(self, "Processing event %i..." % event.id)
    event.exon_bounds = 0
    # iterate through the event's members
    for member in event.members: #{
      self.ProcessMember(member)
      if (member.meta_fields['exon_bounds'] > event.exon_bounds): #{
        event.exon_bounds = member.meta_fields['exon_bounds']
      #} end if
    #} end for
    # write the output
    self.output_file.Write(event.FullDataString())
  #} end def

  def ProcessMember(self, member): #{
    DebugMsg(self, "Processing member %s..." % member.IDString())
    member.meta_fields['exon_bounds'] = 0
    # if the member is an internal gap
    if (member.gap and member.GapIsInternal()): #{
      # test the left (downstream) side of the gap
      breakpoint_coord = min(member.align_info_B.genome_start,
        member.align_info_B.genome_end)
      if (self.BreakpointMatchesBoundary(
          member.align_info_B.chrom, 'left', breakpoint_coord)):
        member.meta_fields['exon_bounds'] += 1
      #} end if
      # test the right (upstream) side of the gap
      breakpoint_coord = max(member.align_info_B.genome_start,
        member.align_info_B.genome_end)
      if (self.BreakpointMatchesBoundary(
          member.align_info_B.chrom, 'right', breakpoint_coord)):
        member.meta_fields['exon_bounds'] += 1
      #} end if
    else:
      # if treating an edge gap event as a split event,
      # alignments might need to be reordered to be in contig order
      if (member.gap and
          member.align_info_A.ctg_start > member.align_info_B.ctg_start):
        DebugMsg(self, "Reordering alignments in gapped split event...")
        align_info_A = member.align_info_B
        align_info_B = member.align_info_A
      else:
        align_info_A = member.align_info_A
        align_info_B = member.align_info_B
      #} end if
      # test the first breakpoint coordinate
      if ("+" == align_info_A.Strand()): #{
        breakpoint_side = 'right'
      else:
        breakpoint_side = 'left'
      #} end if
      breakpoint_coord = align_info_A.genome_end
      if (self.BreakpointMatchesBoundary(
          align_info_A.chrom, breakpoint_side, breakpoint_coord)):
        member.meta_fields['exon_bounds'] += 1
      #} end if
      # test the second breakpoint coordinate
      if ("+" == align_info_B.Strand()): #{
        breakpoint_side = 'left'
      else:
        breakpoint_side = 'right'
      #} end if
      breakpoint_coord = align_info_B.genome_start
      if (self.BreakpointMatchesBoundary(
          align_info_B.chrom, breakpoint_side, breakpoint_coord)):
        member.meta_fields['exon_bounds'] += 1
      #} end if
    #} end if
  #} end def

  def BreakpointMatchesBoundary(self, chrom, j_side, j_coord): #{
    if (chrom not in self.exon_bound_coords): #{
      LogMsg(self, "Warning: no exons on chromosome %s" % chrom)
      return False
    #} end if
    DebugMsg(self, "Breakpoint side: %s Breakpoint coord: %i" %
      (j_side, j_coord))
    offset_range = range(-self.options.buffer_size, self.options.buffer_size+1)
    for offset in offset_range: #{
      test_coord = j_coord + offset
      ExtremeDebugMsg(self, "Offset: %i Test coord: %i" % (offset, test_coord))
      if (test_coord in self.exon_bound_coords[chrom][j_side]): #{
        DebugMsg(self, "  Found exon boundary at %i" % test_coord)
        return True
      #} end if
    #} end if
    return False
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class ExonBoundCounterError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Go through the group members in the input file and "
    "for each member check how many of its breakpoint coordinates match up "
    "with exon boundaries.")
  args = [ "LIB", "BARNACLE_FILE", "GENES_FILE", ]
  usage_string = ("%prog " + " ".join(args) + " [ OPTIONS ]")
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-b", "--buffer-size",
                    type="int", metavar="N",
                    help="Count the breakpoint coordinate as matching an exon "
                         "coordinate if it is within Nbp to either side of "
                         "it. [default:%default]")
  # add option for additional annotations files?
  misc_group = OptionGroup(parser, "Miscellaneous Options")
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
  parser.set_defaults(buffer_size=4,
                      force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  input_dir = os.path.dirname(options.barnacle_path)
  CheckFilePath(options.genes_path, "gene annotations", path_errors)
  options.output_dir = GetOutDir(input_dir, "with_exon_bounds")
  CheckDirPath(options.output_dir, "output", path_errors, create=True)
  # set the log-file name
  options.log_file_name = GetLogPath(options.barnacle_path,
    "exon_bounds", options.output_dir)
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
    options.genes_path    = EnsureAbsPath(args[2])
    if (CheckPaths(options)): #{
      try:
        exon_bound_counter = ExonBoundCounterCls(options)
        WriteCommand(exon_bound_counter, sys.argv)
        exon_bound_counter.CountExonBoundaryBreakpoints()
      except (MyError), e:
        ErrMsg("ERROR while counting exon boundary matching breakpoint "
               "coordinates for Barnacle event members:\n  %s" % e)
        return ES_RUN_ERR
      # end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle data file for that library (BARNACLE_FILE); and the path "
      "to a GTF or UCSC gene annotations file (GENES_FILE).")
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
