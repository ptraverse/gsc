#! /usr/bin/env python
"""
repeats.py

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
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
                           RunOverlapCode, ShouldChromUseChr)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, GetFilePath, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls
from common.candidate_group  import CandidateGroupError

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "REPEAT FLAG SUCCESS"
MSG_FAIL = "REPEAT FLAG FAIL"

class RepeatFlaggerCls: #{
  def __init__(self, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    self.output_files = dict()
    self.num_events = 0
    self.num_passing = 0
    self.rna_regex = re.compile(r"_\(r\|sc\|sn\|srp\|t\)?RNA_[0-9]+$")
  #} end def

  def __del__(self): #{
    for output_file in self.output_files.itervalues(): #{
      output_file.Close()
    #} end for
    for repeats_file in self.options.repeat_files: #{
      repeats_file.close()
    #} end for
    CloseLogFile(self)
  #} end def

  def FlagRepeats(self): #{
    LogMsg(self, "Flagging repeat regions in %s Barnacle events..." %
      self.options.lib)
    start_time = time.time()
    # setup output files
    self.SetupFiles()
    # create event coordinates files
    self.CreateEventCoordsFiles()
    for repeats_file in self.options.repeat_files: #{
      # run overlap code
      self.GetRepeatOverlaps(repeats_file)
    #} end for
    try: #{
      # parse overlap code output
      self.ParseOverlapResults()
    except CandidateGroupError, e:
      raise RepeatFlaggerError("error parsing overlap file: %s" % e)
    #} end try
    LogMsg(self, "Total time spent flagging repeats: %s" %
      TimeSpent(start_time))
  #} end def

  def SetupFiles(self): #{
    file_name = os.path.basename(self.options.barnacle_path)
    out_path = os.path.join(self.options.output_dir, file_name)
    if (self.options.filter_repeats): #{
      out_path = out_path.replace(".pass", ".data")
      if (not out_path.endswith(".data")): #{
        out_path += ".data"
      #} end if
    #} end if
    DebugMsg(self, "Writing events to %s\n" % out_path)
    fail_msg = "cannot create events output file"
    self.output_files["all events"] = FileBoxCls(out_path, "w", fail_msg)
    if (self.options.filter_repeats): #{
      out_path_pass = out_path.replace(".data", ".pass")
      DebugMsg(self, "Writing passing events to %s\n" % out_path_pass)
      fail_msg = "cannot create events output file"
      self.output_files["passing"] = FileBoxCls(out_path_pass, "w", fail_msg)
    #} end if
    for repeats_file in self.options.repeat_files: #{
      repeats_file.log_info = self.log_info
    #} end for
  #} end def

  def CreateEventCoordsFiles(self): #{
    if (self.options.coords_exist or self.options.overlaps_exist): #{
      return
    #} end if
    start_time = time.time()
    LogMsg(self, "Creating event coordinates files...")
    # open the coordinates file
    fail_msg = "cannot create event coordinates file"
    event_coords_file = FileBoxCls(self.options.event_coords_path, "w",
      fail_msg)
    chr_event_coords_file = FileBoxCls(self.options.chr_event_coords_path,
      "w", fail_msg)
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    for event in parser: #{
      DebugMsg(self, "EVENT: %s" % event.DataString())
      for member in event.members: #{
        DebugMsg(self, "MEMBER: %s" % member.DataString())
        # write the coordinates for the event
        self.WriteEventCoords(member, event_coords_file, chr_event_coords_file)
      #} end for
    #} end for
    parser.CloseDataFile()
    event_coords_file.close()
    chr_event_coords_file.close()
    LogMsg(self, "Time spent creating event coordinates files: %s" %
      TimeSpent(start_time))
    if (0 == os.path.getsize(self.options.event_coords_path) or
        0 == os.path.getsize(self.options.chr_event_coords_path)):
      raise RepeatFlaggerError("Event coordinates file is empty!!")
    #} end if
  #} end def

  def WriteEventCoords(self, member, coords_file, chr_coords_file): #{
    DebugMsg(self, "Alignment topology: %s" % member.topology)
    if (member.gap and member.GapIsInternal()): #{
      member_coords = self.GetGapEventCoords(member)
      coords_list = member_coords.CoordsList(use_chr=False)
      chr_coords_list = member_coords.CoordsList(use_chr=True)
    else:
      (member_coords_A, member_coords_B) = self.GetSplitEventCoords(member)
      coords_list = member_coords_A.CoordsList(use_chr=False)
      coords_list.extend(member_coords_B.CoordsList(use_chr=False))
      chr_coords_list = member_coords_A.CoordsList(use_chr=True)
      chr_coords_list.extend(member_coords_B.CoordsList(use_chr=True))
    #} end if
    coords_file.WriteLine("%s" % "\n".join(coords_list))
    chr_coords_file.WriteLine("%s" % "\n".join(chr_coords_list))
  #} end def

  def GetGapEventCoords(self, member): #{
    DebugMsg(self, "Getting gap member coordinates...")
    member_coords = RepeatCoordsCls(member.IDString(), "A",
      log_info=self.log_info)
    # get gap member coordinates
    member_coords.GetGapCoords(member.align_info_B, member.blocks_A,
      member.blocks_B, self.options.region_size)
    DebugMsg(self, "Gap Region: %s\nMerging coords..." %
      (", ".join(member_coords.CoordsList())))
    member_coords.MergeGapCoords()
    DebugMsg(self, "Gap Region: %s" % (", ".join(member_coords.CoordsList())))
    return member_coords
  #} end def

  def GetSplitEventCoords(self, member): #{
    DebugMsg(self, "Getting split member coordinates...")
    region_size = self.options.region_size
    # if treating a gap event as a split event, there is no contig overlap
    if (not member.gap and
        "N/A" != member.meta_fields['ctg_overlap'] and
        0 < member.meta_fields['ctg_overlap']):
      region_size += member.meta_fields['ctg_overlap']
    #} end if
    DebugMsg(self, "Region Size: %i" % region_size)
    # if treating a gap event as a split event,
    # alignments might need to be reordered to be in contig order
    if (member.gap and
        member.align_info_A.ctg_start > member.align_info_B.ctg_start):
      DebugMsg(self, "Reordering alignments in gapped split event...")
      align_info_A = member.align_info_B
      blocks_A = member.blocks_B
      align_info_B = member.align_info_A
      blocks_B = member.blocks_A
    else:
      align_info_A = member.align_info_A
      blocks_A = member.blocks_A
      align_info_B = member.align_info_B
      blocks_B = member.blocks_B
    #} end if
    # get region A member coordinates
    member_coords_A = RepeatCoordsCls(member.IDString(), "A",
      log_info=self.log_info)
    if (align_info_A.genome_start < align_info_A.genome_end): #{
      breakpoint_side = "right"
    else:
      breakpoint_side = "left"
    #} end if
    member_coords_A.GetSplitCoords(align_info_A.chrom, blocks_A,
      breakpoint_side, region_size)
    # get region B member coordinates
    member_coords_B = RepeatCoordsCls(member.IDString(), "B",
      log_info=self.log_info)
    if (align_info_B.genome_start < align_info_B.genome_end): #{
      breakpoint_side = "left"
    else:
      breakpoint_side = "right"
    #} end if
    member_coords_B.GetSplitCoords(align_info_B.chrom, blocks_B,
      breakpoint_side, region_size)
    DebugMsg(self, "Region A: %s\nRegion B: %s" %
      (", ".join(member_coords_A.CoordsList()),
       ", ".join(member_coords_B.CoordsList())))
    return (member_coords_A, member_coords_B)
  #} end def

  def GetRepeatOverlaps(self, repeats_file): #{
    if (self.options.overlaps_exist): #{
      return
    #} end if
    LogMsg(self, "Running overlap code for %s..." % repeats_file.in_path)
    start_time = time.time()
    repeats_path = repeats_file.in_path
    use_chr = ShouldChromUseChr(1, repeats_file.in_path,
      "repeats coordinates", self.log_info)
    if (use_chr): #{
      event_coords_path = self.options.chr_event_coords_path
    else:
      event_coords_path = self.options.event_coords_path
    #} end if
    overlaps_path = repeats_file.OutPath()
    RunOverlapCode(repeats_path, event_coords_path, overlaps_path,
      dpt=self.options.dpt)
    LogMsg(self, "Time spent running overlap code: %s" % TimeSpent(start_time))
  #} end def

  def ParseOverlapResults(self): #{
    start_time = time.time()
    LogMsg(self, "Parsing overlaps files...")
    # create event parser
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    # for each event, add repeats from each repeats file
    for event in parser: #{
      self.AddRepeatsToEvent(event)
      if (self.options.check_stRNA): #{
        self.CheckStructuralRNA(event)
      #} end if
      self.OutputEvent(event)
      #ExtremeDebugMsg(self, "  RNA: %s" % event.RNAString())
    #} end for
    parser.CloseDataFile()
    # close all the repeats files once done with them
    for repeats_file in self.options.repeat_files: #{
      repeats_file.close()
    #} end for
    LogMsg(self, "Time spent parsing overlaps files: %s" %
      TimeSpent(start_time))
  #} end def

  def AddRepeatsToEvent(self, event): #{
    DebugMsg(self, "Adding repeats to event: %i" % event.id)
    #} end if
    if (not self.options.keep_current): #{
      event.ClearRepeats()
    #} end if
    for repeats_file in self.options.repeat_files: #{
      self.AddRepeatsFromFile(event, repeats_file)
    #} end for
  #} end def

  def AddRepeatsFromFile(self, event, repeats_file): #{
    DebugMsg(self, "Getting repeats from %s" % repeats_file.OutPath())
    repeats_dict = repeats_file.GetRepeats(event.id)
    # if the current event overlaps repeats
    if (0 < len(repeats_dict)): #{
      event.AddRepeatsToGroup(repeats_dict, keep_current=True)
    else:
      DebugMsg(self, "Event overlaps no repeats from current file")
    #} end if
  #} end def

  def CheckStructuralRNA(self, event): #{
    #ExtremeDebugMsg(self, "Checking structural RNA for event")
    event.rna = False
    #ExtremeDebugMsg(self, "Setting event.rna to False")
    for member in event.members: #{
      for repeat in (member.repeats_A + member.repeats_B): #{
        if (None != self.rna_regex.search(repeat)): #{
          event.rna = True
          #ExtremeDebugMsg(self, "Setting event.rna to True")
          return
        #} end if
      #} end for
    #} end for
  #} end def

  def OutputEvent(self, event): #{
    if (self.options.filter_repeats and event.OverlapsRepeat()): #{
      event.fail_reasons.add("REPEAT_SEQ")
    #} end if
    event_string = event.FullDataString()
    #ExtremeDebugMsg(self, "Event ouput string: %s" % event_string)
    #ExtremeDebugMsg(self, "  RNA: %s" % event.RNAString())
    # write the event to the all events file
    self.output_files['all events'].Write(event_string)
    # if the event passed the filters
    if (self.options.filter_repeats and event.PassedFilters()): #{
      # write the event to the passing events file
      self.output_files['passing'].Write(event_string)
    #} end if
  #} end def
#} end class

class RepeatCoordsCls: #{
  def __init__(self, full_member_id, region_id, log_info=None): #{
    self.log_info = log_info
    self.full_member_id = full_member_id
    self.region_id = region_id
    self.coords = list()
  #} end def

  def GetGapCoords(self, realignment_coords, original_blocks,
      realignment_blocks, region_size):
    self.chrom = realignment_coords.chrom
    # order the coords
    left  = min(realignment_coords.genome_start, realignment_coords.genome_end)
    right = max(realignment_coords.genome_start, realignment_coords.genome_end)
    # check that the gap-coords are in order
    #if (coords.left > coords.right): #{
    #  raise RepeatFlaggerError("Gap coordinates are not in order: "
    #    "left=%i, right=%i" % (coords.left, coords.right))
    #} end if
    # add the gap coordinates
    #DebugMsg(self, "Initializing coordinates to %i-%i" % (left, right))
    #self.coords.append((left, right))
    #gap_size = coords.Span()
    #gap_size = (coords.right - coords.left) + 1
    DebugMsg(self, "Adding initial blocks from gap region...")
    extra_bases = self.AddCoordsFromBlocks(realignment_blocks,
      realignment_coords.GenomeSpan(), "none")
    bases_used = realignment_coords.GenomeSpan() - extra_bases
    region_size -= bases_used
    #region_size -= gap_size
    region_size = region_size / 2
    DebugMsg(self, "Initial blocks used %ibp, %i bases remaining" %
      (bases_used, region_size))
    if (5 > region_size): #{
      return
    #} end if
    # find the last block before the gap
    before_gap = list()
    after_gap = list()
    for block in original_blocks: #{
      (block_left, block_right) = block
      if (block_left > block_right): #{
        raise RepeatFlaggerError("Block coordinates are not in order: "
          "left=%i, right=%i" % (block_left, block_right))
      #} end if
      #(block_left, block_right) = map(int, block.split("-"))
      if (block_right < left): #{
        before_gap.append((block_left, block_right))
        #before_gap.append("%i-%i" % (block_left, block_right))
      elif (block_left < left): #{
        before_gap.append((block_left, left))
        #before_gap.append("%i-%i" % (block_left, realignment_coords.left))
      #} end if
      if (block_left > right): #{
        after_gap.append((block_left, block_right))
        #after_gap.append("%i-%i" % (block_left, block_right))
      elif (block_right > right): #{
        after_gap.append((right, block_right))
        #after_gap.append("%i-%i" % (realignment_coords.right, block_right))
      #} end if
    #} end for
    DebugMsg(self, "Adding blocks before gap region...")
    self.AddCoordsFromBlocks(reversed(before_gap), region_size, "right")
    DebugMsg(self, "Adding blocks after gap region...")
    self.AddCoordsFromBlocks(after_gap, region_size, "left")
  #} end def

  def GetSplitCoords(self, chrom, blocks, breakpoint_side, region_size): #{
    self.chrom = chrom
    if ("left" == breakpoint_side): #{
      # start from the left, and add block coords until
      # the desired number of bases is reached
      self.AddCoordsFromBlocks(blocks, region_size, "left")
    elif ("right" == breakpoint_side): #{
      # start from the right, and add block coords until
      # the desired number of bases is reached
      self.AddCoordsFromBlocks(reversed(blocks), region_size, "right")
    else:
      raise RepeatFlaggerError("invalid breakpoint side: %s" % breakpoint_side)
    #} end if
  #} end def

  def AddCoordsFromBlocks(self, blocks, num_bases, from_side): #{
    DebugMsg(self, "Adding coordinates from blocks, "
      "starting on the %s..." % from_side)
    for block in blocks: #{
      if (0 >= num_bases): #{
        break
      #} end if
      (block_left, block_right) = block
      #(block_left, block_right) = map(int, block.split("-"))
      block_len = (block_right - block_left) + 1
      #print "Block Len: %i, Bases Remaining: %i" % (block_len, num_bases)
      if (block_len > num_bases): #{
        if("left" == from_side): #{
          # add just the left side of the block
          block_right = (block_left + num_bases) - 1
        elif("right" == from_side): #{
          # add just the right side of the block
          block_left = (block_right - num_bases) + 1
        else:
          raise RepeatFlaggerError("invalid from side: %s" % from_side)
        #} end if
        block_len = (block_right - block_left) + 1
      #} end if
      # add the block
      DebugMsg(self, "Adding block: %i-%i" % (block_left, block_right))
      self.coords.append((block_left, block_right))
      # reduce the number of bases remaining
      num_bases -= block_len
    #} end for
    return num_bases
  #} end def

  def MergeGapCoords(self): #{
    DebugMsg(self, "Merging gap coordinates...")
    # if there is only a single set of coordinates in the list,
    # no merging is necessary
    if (2 > len(self.coords)): #{
      return
    #} end if
    # sort the coordinates by their left side
    self.coords.sort(key=lambda block: block[0])
    new_coords = list()
    (curr_left, curr_right) = self.coords[0]
    for coords in self.coords[1:]: #{
      (next_left, next_right) = coords
      DebugMsg(self, "  CURR: %i-%i\n" % (curr_left, curr_right) +
        "  NEXT: %i-%i" % (next_left, next_right))
      if (next_left > curr_right + 1): #{
        new_coords.append((curr_left, curr_right))
        (curr_left, curr_right) = (next_left, next_right)
      else:
        curr_right = max(curr_right, next_right)
      #} end if
    #} end for
    new_coords.append((curr_left, curr_right))
    self.coords = new_coords
  #} end def

  def CoordsList(self, use_chr=True): #{
    return [self.ToString(i, use_chr) for i in range(len(self.coords))]
  #} end def

  def ToString(self, index, use_chr=True): #{
    if (use_chr): #{
      chrom = "chr%s" % self.chrom
    else:
      chrom = self.chrom
    #} end if
    data_list = list([
      chrom,
      "%i" % self.coords[index][0],
      "%i" % self.coords[index][1],
      "%s%s" % (self.full_member_id, self.region_id),
    ])
    return " ".join(data_list)
  #} end def
#} end class

class RepeatCoordsFileCls: #{
  def __init__(self, in_path): #{
    self.in_path  = EnsureAbsPath(in_path)
    self.out_file = None
    self.event_id = -1
  #} end def

  def __del__(self): #{
    self.close()
  #} end def

  def Setup(self, lib_name, out_dir): #{
    in_file_name  = os.path.splitext(os.path.basename(self.in_path))[0]
    DebugMsg(self, "In file name: %s" % in_file_name)
    out_file_name = "%s.%s.olap" % (lib_name, in_file_name)
    DebugMsg(self, "Out file name: %s" % out_file_name)
    out_path      = os.path.join(out_dir, out_file_name)
    fail_msg      = "cannot read overlaps file"
    self.out_file = FileBoxCls(out_path, 'r', fail_msg)
  #} end def

  def OutPath(self): #{
    return self.out_file.path
  #} end def

  def GetRepeats(self, event_id): #{
    repeats_dict = dict()
    try: #{
      # ignore groups before the current event
      while (self.event_id < event_id): #{
        self.ParseOverlapLine()
      #} end while
      # process all overlap lines for the current event
      while (self.event_id == event_id): #{
        DebugMsg(self, "Adding repeats: %s" % self.repeats)
        # ensure the repeats dictionary is set up
        if (self.member_id not in repeats_dict): #{
          repeats_dict[self.member_id] = dict()
        #} end if
        repeats_dict_for_member = repeats_dict[self.member_id]
        if (self.region_id not in repeats_dict_for_member): #{
          repeats_dict_for_member[self.region_id] = dict()
        #} end if
        repeats_dict_for_region = repeats_dict_for_member[self.region_id]
        for repeat_flag in self.GetRepeatFlags(self.repeats): #{
          repeats_dict_for_region[repeat_flag.lower()] = repeat_flag
        #} end for
        self.ParseOverlapLine()
      #} end while
    except StopIteration,e:
      DebugMsg(self, "End of overlaps file reached!")
    #} end try
    return repeats_dict
  #} end def

  def ParseOverlapLine(self): #{
    overlap_line = self.out_file.next()
    DebugMsg(self, "OVERLAP: %s" % overlap_line)
    (chrom, left, right, id_list, repeats) = overlap_line.split(" ", 4)
    id_match = re.search(
      r"(?P<event>\d+)(?P<member>[a-z]+)\)?(?P<region>[AB])", id_list)
    if (None == id_match): #{
      raise RepeatFlaggerError("could not get event and region IDs from "
         "overlap line: \"%s\"\n  %s" % (id_list, overlap_line))
    #} end if
    DebugMsg(self, "GROUP ID: %s, MEMBER: %s, REGION: %s" %
      (id_match.group('event'), id_match.group('member'),
       id_match.group('region')))
    #} end if
    self.event_id = int(id_match.group('event'))
    self.member_id = id_match.group('member')
    self.region_id = id_match.group('region')
    self.repeats = repeats
  #} end def

  def GetRepeatFlags(self, repeats_str): #{
    repeat_flags = list()
    # remove any leading "^" and then split on remaining "^"s
    repeats = repeats_str.strip("^").split("^")
    # get the repeat flags for each repeat found
    for repeat_str in repeats: #{
      # remove the XcountX information
      repeat_str = re.sub(r"X\d+X$", "", repeat_str)
      # split the id parts
      #id_cols = repeat_str.split(" ")
      #flag_parts = list()
      #for col in self.cols_to_use: #{
      #  flag_parts.append(id_cols[col])
      ##} end for
      #repeat_flag = "/".join(flag_parts)
      #repeat_flags.append(repeat_flag)
      repeat_flags.append(repeat_str.replace(" ", "/"))
    #} end for
    return repeat_flags
  #} end def

  def close(self): #{
    if (None != self.out_file): #{
      self.out_file.close()
    #} end if
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class RepeatFlaggerError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Flags events when the event-region overlaps an "
    "annotated repeat sequence")
  args = [ "LIB", "BARNACLE_FILE", "REPEAT_FILES", ]
  usage_string = "\n".join(["%prog " + " ".join(args) + " [ OPTIONS ]",
    "    REPEAT_FILES can be a single path, or a comma-delimited list."])
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--region-size",
                    type="int", metavar="N",
                    help="Look for repeats in a region N bp wide around the "
                         "event coordinates. [default: %default]")
  parser.add_option("--repeat-filter-off",
                    action="store_false", dest="filter_repeats",
                    help="Do not filter out events that are flagged as "
                         "overlapping repeat sequence. [default]")
  parser.add_option("--repeat-filter-on",
                    action="store_true", dest="filter_repeats",
                    help="Filter out events that are flagged as "
                         "overlapping repeat sequence.")
  parser.add_option("--output-dir",
                    metavar="DIR",
                    help="Put output files in directory DIR.")
  parser.add_option("--coords-exist",
                    action="store_true",
                    help="Use an already existing coordinates file.")
  parser.add_option("--overlaps-exist",
                    action="store_true",
                    help="Use an already existing overlaps file.")
  parser.add_option("--replace-current-repeats",
                    action="store_false", dest="keep_current",
                    help="Replace any repeats the events are already "
                         "annotated with with the new repeats. [default]")
  parser.add_option("--keep-current-repeats",
                    action="store_true", dest="keep_current",
                    help="Add the new repeats to any repeats the events "
                         "are already annotated with.")
  parser.add_option("--rna", "--structural-RNA-check",
                    action="store_true", dest="check_stRNA",
                    help="Check whether the event overlaps a structural RNA. "
                      "[default]")
  parser.add_option("--no-rna", "--no-structural-RNA-check",
                    action="store_false", dest="check_stRNA",
                    help="Do not check whether the event overlaps a "
                      "structural RNA.")
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  misc_group.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(region_size=100,
                      filter_repeats=False,
                      coords_exist=False,
                      overlaps_exist=False,
                      keep_current=False,
                      check_stRNA=True,
                      dpt=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  # get the output directory
  options.output_dir = GetOutDir(os.path.dirname(options.barnacle_path),
    "with_repeats")
  options.event_coords_path = GetFilePath(options.output_dir,
    options.lib, "repeats.coords")
  if (options.coords_exist): #{
    if (not os.path.exists(options.event_coords_path) and
        not options.overlaps_exist):
      path_errors.append("Cannot find event coordinates file, do not use "
        "\"--coords-exist\" option and a new event coordinates file will "
        "be created: %s" % options.event_coords_path)
    #} end if
  else:
    if (os.path.exists(options.event_coords_path)): #{
      path_errors.append("Event coordinates file already exists, either "
        "remove it or use \"--coords-exist\" option to use existing event "
        "coordinates file: %s" % options.event_coords_path)
    #} end if
  #} end if
  options.chr_event_coords_path = GetFilePath(options.output_dir,
    options.lib, "repeats.chr.coords")
  if (options.coords_exist): #{
    if (not os.path.exists(options.chr_event_coords_path) and
        not options.overlaps_exist):
      path_errors.append("Cannot find event coordinates file with \"chr\", "
        "do not use \"--coords-exist\" option and a new event coordinates "
        "file will be created: %s" % options.chr_event_coords_path)
    #} end if
  else:
    if (os.path.exists(options.chr_event_coords_path)): #{
      path_errors.append("Event coordinates file with \"chr\" already exists, "
        "either remove it or use \"--coords-exist\" option to use existing "
        "event coordinates file: %s" % options.chr_event_coords_path)
    #} end if
  #} end if
  for repeats_file in options.repeat_files: #{
    CheckFilePath(repeats_file.in_path, "repeat coordinates", path_errors)
    repeats_file.Setup(options.lib, options.output_dir)
    if (options.overlaps_exist): #{
      if (not os.path.exists(repeats_file.OutPath())): #{
        path_errors.append("Cannot find overlaps file, do not use "
          "\"--overlaps-exist\" option and a new overlaps file will be "
          "created: %s" % repeats_file.OutPath())
      #} end if
    else:
      if (os.path.exists(repeats_file.OutPath())): #{
        path_errors.append("Overlaps file already exists, either remove it "
          "or use \"--overlaps-exist\" option to use existing overlaps "
          "file: %s" % repeats_file.OutPath())
      #} end if
    #} end if
  #} end for
  # check the output directory
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors, create=True)
    # get the log-file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "repeats", options.output_dir)
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
    options.repeat_files  = [RepeatCoordsFileCls(path) for
      path in args[2].strip(",").split(",")]
    if (CheckPaths(options)): #{
      try: #{
        ac_repeat_flagger = RepeatFlaggerCls(options)
        WriteCommand(ac_repeat_flagger, sys.argv)
        ac_repeat_flagger.FlagRepeats()
      except (MyError), e:
        ErrMsg("ERROR while flagging repeat sequences in BARANCLE events:"
          "\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle data file for that library (BARNACLE_FILE); and a "
      "comma-delimited list of paths to BED files containing repeat "
      "element coordinates (REPEAT_FILES).")
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
