#! /usr/bin/env python
"""
predict_events.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# TODO:
# Output lib_count in "across all libraries" genes list
# implement --filter-members option
#   (add this to quick_filter and bubble_up too)
# merge "annotations" and "filter-gene-names" options

# DONE:

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
  ShouldChromUseChr, RunOverlapCode, RemoveIndexingFlag, CheckConfigCommands)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls
from fusion_predictor import FusionPredictorCls
from ptd_predictor import PTDPredictorCls
from itd_predictor import ITDPredictorCls
from realigner import RealignerCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "EVENT PREDICTION SUCCESS"
MSG_FAIL = "EVENT PREDICTION FAIL"

class EventPredictionCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    if (not hasattr(self.options, "realign")): #{
      self.options.realign = False
    #} end if
    if (self.options.realign): #{
      CheckConfigCommands(self, "blat")
    #} end if
    self.predictors = dict()
    if (self.options.predict_fusions): #{
      predictor = FusionPredictorCls(options, log_info=self.log_info)
      self.predictors[predictor.key] = predictor
    #} end if
    if (self.options.predict_ptds): #{
      predictor = PTDPredictorCls(options, log_info=self.log_info)
      self.predictors[predictor.key] = predictor
    #} end if
    if (self.options.predict_itds): #{
      predictor = ITDPredictorCls(options, log_info=self.log_info)
      self.predictors[predictor.key] = predictor
    #} end if
    self.use_chr = False
    #self.postpone_gene_check = False
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def PredictEvents(self): #{
    LogMsg(self, "Predicting events...")
    start = time.time()
    # get the reference gene names, if a path is given
    #self.ref_gene_names = GetGeneNamesFromFile(self.options.gene_names_path,
    #  self.log_info)
    group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    # recheck breakpoint exons
    self.RecheckBreakpointExons(group_parser)
    realigner = None
    if (self.options.realign): #{
      realigner = RealignerCls(self.options, self.log_info)
    #} end if
    # potential_events[bio_type][group_id] = event and gene sets object
    #potential_events = dict([(predictor.key, dict()) for
    #  predictor in self.predictors])
    LogMsg(self, "Processing candidate groups...")
    process_start = time.time()
    for group in group_parser: #{
      # get the breakpoint exons for the group
      self.GetBreakpointExons(group)
      # check whether the event is any biologically typed event
      #self.CheckEvent(group, output_files, lib_info.lib_name, potential_events)
      # attempt to predict events of each specified type
      # from the current candidate group
      for predictor in self.predictors.itervalues(): #{
        good_members = list()
        if (predictor.ProcessGroup(group, good_members) and
            None != realigner): #{
          #realigner.UpdateContigs(group, good_members, predictor.store_seq)
          realigner.UpdateContigs(group, good_members, predictor.key)
        #} end if
      #} end for
    #} end for
    LogMsg(self, "Time spent processing candidate groups: %s" %
      TimeSpent(process_start))
    if ("itd" in self.predictors and
        0 < self.predictors["itd"].num_over_aligned): #{
      LogMsg(self, "WARNING: %i gap candidates have aligned length greater "
        "than gap length!" % self.predictors["itd"].num_over_aligned)
    #} end if
    #if ('event_coords' in output_files): #{
    #  output_files['event_coords'].Close()
    #  self.RecheckExonOverlap(output_files, potential_events, lib_info.lib_name)
    #} end if
    if (None != realigner and 0 < len(realigner.contigs)): #{
      realigner.RealignContigs()
      LogMsg(self, "Before realignment:")
      for predictor in self.predictors.itervalues(): #{
        LogMsg(self, "  Number of %s predictions: %i" %
          (predictor.description, predictor.num_predictions))
        if (0 == predictor.num_predictions): #{
          continue
        #} end if
        if ("itd" in predictor.key or "fusion" in predictor.key): #{
          predictor.LoadTranscriptSequences(realigner.contigs)
        #} end if
        predictor.ReprocessPredictions(realigner.contigs)
        #predictor.ReprocessPredictions(realigner.contigs,
        #  realigner.contig_seqs)
      #} end for
      LogMsg(self, "%s\nAfter realignment:" % ("-"*40))
    #} end if
    for predictor in self.predictors.itervalues(): #{
      LogMsg(self, "Number of %s predictions: %i" %
        (predictor.description, predictor.num_predictions))
    #} end for
    LogMsg(self, "Time spent predicting events: %s" % TimeSpent(start))
  #} end def

  #def CreateOutputFiles(self, input_path): #{
  #  input_file_name = os.path.basename(input_path)
  #  input_root = os.path.splitext(input_file_name)[0]
  #  output_files = dict()
  #  # setup the coordinates file for rechecking exon overlaps
  #  self.SetupEventCoordsFile(input_root, output_files)
  #  return output_files
  #} end def

  def RecheckBreakpointExons(self, group_parser): #{
    if (None == self.options.breakpoint_exons): #{
      self.overlaps_file = None
      return
    #} end if
    if (self.options.use_existing_group_coords): #{
      LogMsg(self, "Using existing group coordinates file.")
    else:
      group_coords_file = self.CreateGroupCoordsFile()
      for group in group_parser: #{
        #ExtremeDebugMsg(self, "Writing coordinates for group %i" % group.id)
        self.WriteGroupCoords(group, group_coords_file)
      #} end for
      group_parser.Close()
      group_coords_file.Close()
    #} end if
    self.CreateOverlapsFile()
  #} end def

  #def SetupEventCoordsFile(self, input_root, output_files): #{
  def CreateGroupCoordsFile(self): #{
    # check whether to use "chr" in chromosome names in coordinates file
    self.use_chr = ShouldChromUseChr(1, self.options.breakpoint_exons,
      "exon coordinates", self.log_info)
    # open the group coordinates file
    #output_files['event_coords'] = FileBoxCls(group_coords_path, "w",
    group_coords_file = FileBoxCls(self.options.group_coords_path, "w",
      "cannot create event coordinates file")
    #self.postpone_gene_check = True
    return group_coords_file
  #} end def

  #def WriteEventCoords(self, event, group_coords_file): #{
  def WriteGroupCoords(self, event, group_coords_file): #{
    for member in event.members: #{
      if (member.gap): #{
        # write gap event coordinates
        self.WriteGapGroupCoords(member, group_coords_file)
      else:
        # write split event coordinates
        self.WriteSplitGroupCoords(member, group_coords_file)
      #} end if
    #} end for
  #} end def

  def WriteGapGroupCoords(self, member, group_coords_file): #{
    gap_coords = GroupCoordsCls(
      member.align_info_B.chrom,
      min(member.align_info_B.genome_start, member.align_info_B.genome_end),
      max(member.align_info_B.genome_start, member.align_info_B.genome_end),
      "%sA" % member.IDString(),
      self.use_chr
    )
    group_coords_file.WriteLine("%s" % gap_coords.ToString())
  #} end def

  def WriteSplitGroupCoords(self, member, group_coords_file): #{
    split_coords_A = GroupCoordsCls(
      member.align_info_A.chrom,
      member.align_info_A.genome_end - self.options.event_buffer,
      member.align_info_A.genome_end + self.options.event_buffer,
      "%sA" % member.IDString(),
      self.use_chr
    )
    group_coords_file.WriteLine("%s" % split_coords_A.ToString())
    split_coords_B = GroupCoordsCls(
      member.align_info_B.chrom,
      member.align_info_B.genome_start - self.options.event_buffer,
      member.align_info_B.genome_start + self.options.event_buffer,
      "%sB" % member.IDString(),
      self.use_chr
    )
    group_coords_file.WriteLine("%s" % split_coords_B.ToString())
  #} end def

  #def RecheckExonOverlap(self, output_files, potential_events, lib_name): #{
  #  LogMsg(self, "Rechecking exon overlap...")
  #  start = time.time()
  #  # run overlap code
  #  overlaps_path = self.RunOverlapCode(output_files['group_coords'].path)
  #  try: #{
  #    # parse overlap code output
  #    self.ParseOverlapResults(overlaps_path, potential_events)
  #  except ACEventGroupError, e:
  #    raise EventPredictionError("error parsing overlap file: %s" % e)
  #  #} end try
  #  self.ProcessPotentialEvents(potential_events, output_files, lib_name)
  #  LogMsg(self, "Time spent rechecking exon overlaps: %s" % TimeSpent(start))
  #} end def

  def CreateOverlapsFile(self): #{
    if (self.options.use_existing_overlaps): #{
      LogMsg(self, "Using existing breakpoint/transcript overlaps file.")
    else:
      LogMsg(self, "Running overlap code...")
      if (hasattr(self, "log_file") and None != self.log_file): #{
        self.log_file.Flush()
      #} end if
      start = time.time()
      RunOverlapCode(self.options.breakpoint_exons,
        self.options.group_coords_path, self.options.overlaps_path,
        dpt=self.options.dpt)
      LogMsg(self, "Time spent running overlap code: %s" % TimeSpent(start))
    #} end if
    self.overlaps_file = FileBoxCls(self.options.overlaps_path, "r",
      "cannot read exon/group overlaps file")
    #self.GetNextExonOverlap()
  #} end def

  def GetBreakpointExons(self, group): #{
    if (not hasattr(self, "overlaps_file") or None == self.overlaps_file): #{
      return
    #} end if
    ExtremeDebugMsg(self, "Getting breakpoint exons for group %i" % group.id)
    # clear any previous breakpoint genes
    group.ClearBPGenes()
    if (not hasattr(self, "curr_overlap")): #{
      self.curr_overlap = None
    #} end if
    # skip overlaps for groups that come before the current group
    while (None == self.curr_overlap or
        self.curr_overlap.group_id < group.id): #{
      try: #{
        self.GetNextExonOverlap()
      except StopIteration:
        return
      #} end try
    #} end while
    # create a dictionary of the members of the current group
    members_dict = dict()
    for member in group.members: #{
      members_dict[member.candidate_id] = member
    #} end for
    # get all overlaps for the current group
    while (self.curr_overlap.group_id == group.id): #{
      if (self.curr_overlap.member_id in members_dict): #{
        self.AddBreakPointGene(members_dict[self.curr_overlap.member_id])
      #} end if
      try: #{
        self.GetNextExonOverlap()
      except StopIteration:
        return
      #} end try
    #} end while
  #} end def

  def GetNextExonOverlap(self): #{
    if (not hasattr(self, "overlaps_file") or None == self.overlaps_file): #{
      ExtremeDebugMsg(self, "Setting current overlap to \"None\".")
      self.curr_overlap = None
      return
    #} end if
    overlap = ExonOverlapCls(self.overlaps_file.next())
    self.curr_overlap = overlap
    ExtremeDebugMsg(self, "Current overlap = G%i%s r%s exons: %s" %
      (overlap.group_id, overlap.member_id, overlap.region_id,
      ",".join(overlap.exons)))
  #} end def

  def AddBreakPointGene(self, member): #{
    if (None == self.curr_overlap): #{
      return
    #} end if
    if (self.curr_overlap.group_id != member.group_id): #{
      raise EventPredictionError("Group ID: %i does not match overlap ID: %i" %
        (member.group_id, self.curr_overlap.group_id))
    #} end if
    if (self.curr_overlap.member_id != member.candidate_id): #{
      raise EventPredictionError("Candidate ID: %s " % member.candidate_id +
        "does not match overlap ID: %s" % self.curr_overlap.member_id)
    #} end if
    member.AddGenes("breakpoint_%s" % self.curr_overlap.region_id,
      self.curr_overlap.exons)
  #} end def
#} end class

class GroupCoordsCls: #{
  def __init__(self, chrom, left, right, id, use_chr): #{
    self.chrom = chrom
    self.left  = left
    self.right = right
    self.id    = id
    self.use_chr = use_chr
  #} end def

  def ToString(self): #{
    if (self.use_chr): #{
      chromosome = "chr%s" % self.chrom
    else:
      chromosome = self.chrom
    #} end if
    data_list = list([
      chromosome,
      "%i" % self.left,
      "%i" % self.right,
      self.id,
    ])
    return " ".join(data_list)
  #} end def
#} end class

class ExonOverlapCls: #{
  def __init__(self, overlap_line, log_info=None): #{
    self.log_info = log_info
    DebugMsg(self, "OVERLAP: %s" % overlap_line)
    # split the line into fields
    (chrom, left, right, id_str, exon_str) = overlap_line.split(" ", 4)
    # get the group, member, and region IDs
    id_patt = r"^(?P<group>\d+)(?P<member>[a-z]+)\)?(?P<region>[AB])$"
    id_match = re.search(id_patt, id_str)
    if (None == id_match): #{
      raise EventPredictionError("could not get ids from overlap line: "
        "\"%s\"\n" % id_str + "  %s" % overlap_line)
    #} end if
    self.group_id = int(id_match.group('group'))
    self.member_id = id_match.group('member')
    self.region_id = id_match.group('region')
    # get the exons overlapped
    exon_str = exon_str.replace(" ", "_")
    self.exons = set([RemoveIndexingFlag(exon) for
      exon in exon_str.lstrip("^").split("^")])
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class EventPredictionError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Predict events of specific types")
  args = [ "LIB", "BARNACLE_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  fusion_grp = OptionGroup(parser, "Fusion Prediction Options")
  fusion_grp.add_option("--predict-fusions",
                    action="store_true", dest="predict_fusions",
                    help="Predict fusion events. [default]")
  fusion_grp.add_option("--no-fusions",
                    action="store_false", dest="predict_fusions",
                    help="Do not predict fusion events")
  fusion_grp.add_option("--with-gene-directions",
                    action="store_true", dest="use_gene_directions",
                    help="Only predict a fusion event if the transcription "
                         "direction of the genes involved is continuous "
                         "across the contig. [default]")
  fusion_grp.add_option("--without-gene-directions",
                    action="store_false", dest="use_gene_directions",
                    help="Do not consider transcription direction of genes "
                         "involved.")
  fusion_grp.add_option("--transcript-annotations",
                    metavar="PATH",
                    help="PATH is the path to a file containing gene "
                         "transcript annotations. Used for gene directions.")
                         #"and ITD realignment.")
  fusion_grp.add_option("--use-conflicts",
                    action="store_true",
                    help="When a gene has conflicting directions in the "
                         "annotations file, believe that the gene really "
                         "does go both ways. [default]")
  fusion_grp.add_option("--no-use-conflicts",
                    action="store_false", dest="use_conflicts",
                    help="Do not consider genes going in both directions. "
                      "Raise an error if the annotations file contains "
                      "conflicting directions.")
  fusion_grp.add_option("--ignore-conflicts",
                    action="store_true",
                    help="When a gene has conflicting directions in the "
                         "annotations file, just ignore that gene, rather "
                         "than halting the script.")
  fusion_grp.add_option("--include-introns",
                    action="store_true",
                    help="Predict fusions that overlap introns even if they "
                      "do not overlap any coding regions.")
  fusion_grp.add_option("--include-nearby",
                    action="store_true",
                    help="Predict fusions that do not overlap genes if they "
                      "are close to genes.")
  fusion_grp.add_option("--fusion-dup-filter",
                    action="store_true",
                    help="Do not predict a fusion when the contig could be "
                         "explained by a duplication event instead. "
                         "[default]")
  fusion_grp.add_option("--no-fusion-dup-filter",
                    action="store_false", dest="fusion_dup_filter",
                    help="Do not check whether fusion contigs could be "
                         "explained by duplication events instead.")
  fusion_grp.add_option("--min-unique-align",
                    type="int", metavar="N",
                    help="Require that at least Nbp of the contig aligns to "
                      "only one of the fusion partners. [default: %default]")
  fusion_grp.add_option("--min-exon-bounds",
                    type="int", metavar="N",
                    help="Require that N (0, 1, or 2) of the breakpoints "
                      "match up with annotated exon boundaries. [default: "
                      "%default]")
  fusion_grp.add_option("--read-through",
                    type="int", metavar="N",
                    help="Consider colinear split alignments closer than Nbp "
                      "to be read-through events. [default: %default]")
  parser.add_option_group(fusion_grp)
  #parser.add_option("--filter-gene-names",
  #                  dest="gene_names_path", metavar="FILE",
  #                  help="Use only gene names found in FILE (should be "
  #                       "tab-separated file, with gene name in column 13).")
  ptd_group = OptionGroup(parser, "PTD Prediction Options")
  ptd_group.add_option("--predict-PTDs",
                    action="store_true", dest="predict_ptds",
                    help="Predict PTD events. [default]")
  ptd_group.add_option("--no-PTDs",
                    action="store_false", dest="predict_ptds",
                    help="Do not predict PTD events")
  parser.add_option_group(ptd_group)
  itd_group = OptionGroup(parser, "ITD Prediction Options")
  itd_group.add_option("--predict-ITDs",
                    action="store_true", dest="predict_itds",
                    help="Predict ITD events. [default]")
  itd_group.add_option("--no-ITDs",
                    action="store_false", dest="predict_itds",
                    help="Do not predict ITD events")
  itd_group.add_option("--allow-ITD-repeats",
                    action="store_false", dest="filter_itd_repeats",
                    help="Allow ITD event prediction for candidate "
                         "contigs that overlap annotated repeats.")
  itd_group.add_option("--remove-ITD-repeats",
                    action="store_true", dest="filter_itd_repeats",
                    help="Do not predict an ITD event if the candidate "
                         "contig overlaps any annotated repeats. [default]")
  itd_group.add_option("--require-internal-gaps",
                    action="store_true", dest="require_internal_gaps",
                    help="Require that the duplicated sequence is internal "
                         "to the contig for gapped ITD "
                         "event predictions.")
  itd_group.add_option("--allow-edge-gaps",
                    action="store_false", dest="require_internal_gaps",
                    help="Allow the duplicated sequence to be at the very "
                         "edge of the contig for gapped partial exon "
                         "duplication predictions. [default]")
  itd_group.add_option("--min-edge-gap-fraction",
                    type="float", metavar="F",
                    help="If the gap is at the very edge of a contig, "
                         "require that the fraction of the gap involved "
                         "in the duplication is at least F. This option "
                         "is ignored if the --require-internal-gaps option "
                         "is used. [default: %default]")
  itd_group.add_option("--allow-non-gap-ITDs",
                    action="store_true", dest="allow_non_gap_itds",
                    help="Look for ITDs in groups with junction duplication, "
                         "end duplication, and non-colinear topologies "
                         "as well as gap topologies.")
  itd_group.add_option("--only-gap-ITDs",
                    action="store_false", dest="allow_non_gap_itds",
                    help="Only look for ITDs in groups with gap "
                      "topologies [default].")
  itd_group.add_option("--exclude-non-coding",
                    action="store_true",
                    help="Do not predict ITD events in non-coding genes. "
                         "[default]")
  itd_group.add_option("--allow-non-coding",
                    action="store_false", dest="exclude_non_coding",
                    help="Report ITD events predicted in non-coding genes.")
  parser.add_option_group(itd_group)
  parser.add_option("--get-breakpoint-exons",
                    metavar="FILE", dest="breakpoint_exons",
                    help="Use the exon coordinates in FILE to get the exon "
                         "overlapped by the breakpoint in each event (used "
                         "for ITDs).")
  parser.add_option("--event-buffer",
                    type="int", metavar="N",
                    help="For split events, check whether the breakpoint "
                         "coordinates are within N bases of an exon. "
                         "[default: %default]")
  parser.add_option("--use-existing-overlaps",
                    action="store_true",
                    help="If a breakpoint/transcripts overlap file already "
                         "exists, just use it rather than generating a new "
                         "one.")
  parser.add_option("--transcript-sequences",
                    metavar="PATH", dest="tran_seq_path",
                    help="The path to a file containing transcript sequences "
                         "to realign candidate contigs against, for avoiding "
                         "false positives due to problems with contig to "
                         "genome alignments.")
  parser.add_option("--contig-sequences",
                    metavar="PATH", dest="ctg_seq_path",
                    help="The path to a file containing candidate contig "
                         "sequences to realign to transcript sequences, for "
                         "avoiding false positives due to problems with "
                         "contig to genome alignments.")
  parser.add_option("--use-existing-realigns",
                    action="store_true",
                    help="If contig-to-transcript realignment files already "
                         "exist, just use them rather than generating new "
                         "ones.")
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--read-length",
      type="int", metavar="N",
      help="The length of the sequenced reads (used for ITD full duplicate "
        "alignment filter). [default: %default]")
  misc_group.add_option("-p", "--pretty",
                    action="store_true",
                    help="Print more readable output file as well as "
                         "standard data file.")
  misc_group.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  misc_group.add_option("-f", "--force",
                    action="store_true",
                    help="Force prediction to take place, even if the output "
                         "directory already exists.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  misc_group.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(#list_input=True,
                      predict_fusions=True,
                      use_gene_directions=True,
                      use_conflicts=True,
                      ignore_conflicts=False,
                      include_introns=False,
                      include_nearby=False,
                      fusion_dup_filter=True,
                      min_unique_align=5,
                      min_exon_bounds=0,
                      read_through=50000,
                      predict_ptds=True,
                      predict_itds=True,
                      filter_itd_repeats=True,
                      require_internal_gaps=False,
                      min_edge_gap_fraction=0.80,
                      allow_non_gap_itds=False,
                      exclude_non_coding=True,
                      event_buffer=5,
                      use_existing_overlaps=False,
                      use_existing_realigns=False,
                      read_length=75,
                      pretty=False,
                      dpt=False,
                      force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (options.ignore_conflicts): #{
    options.use_conflicts = False
  #} end if
  #if (options.ignore_conflicts and options.use_conflicts): #{
  #  ErrMsg("Cannot use both \"--ignore-conflicts\" and \"--use-conflicts\" "
  #         "options at the same time. Please choose one or the other.")
  #  opts_good = False
  #} end if
  if (0 > options.min_exon_bounds or 2 < options.min_exon_bounds): #{
    ErrMsg("\"--min-exon-bounds\" value must be 0, 1, or 2. Not %i" %
      options.min_exon_bounds)
    opts_good = False
  #} end if
  path_errors = list()
  options.realign = False
  if (None != options.tran_seq_path): #{
    CheckFilePath(options.tran_seq_path, "transcript sequences", path_errors)
    if(None == options.ctg_seq_path): #{
      ErrMsg("You must provide the path to a candidate contig sequences file "
        "when you provide a transcript sequences path.")
      opts_good = False
    #} end if
    options.realign = True
  #} end if
  if (None != options.ctg_seq_path): #{
    CheckFilePath(options.ctg_seq_path, "candidate contig sequences",
      path_errors)
    if(None == options.tran_seq_path): #{
      ErrMsg("You must provide the path to a transcript sequences file "
        "when you provide a candidate contig sequences path.")
      opts_good = False
    #} end if
  #} end if
  CheckFilePath(options.barnacle_path, "groups", path_errors)
  # get the input directory
  (input_dir, file_name) = os.path.split(options.barnacle_path)
  # if checking the gene directions
  if (options.use_gene_directions): #{  or options.realign_itds):
    # ensure a gene annotations path was given
    if (None == options.transcript_annotations): #{
      path_errors.append("You must give the path to a gene transcript "
        "annotations file when using the --with-gene-directions option.")
        #"or --realign-ITDs option.")
    else:
      CheckFilePath(options.transcript_annotations,
        "gene transcript annotations", path_errors)
    #} end if
  #} end if
  if (None != options.breakpoint_exons): #{
    CheckFilePath(options.breakpoint_exons, "exon coordinate annotations",
      path_errors)
  #} end if
  options.output_dir = GetOutDir(input_dir, "predicted_events")
  # get the path for the group coordinates file
  input_file_name = os.path.basename(options.barnacle_path)
  input_root = os.path.splitext(input_file_name)[0]
  options.group_coords_path = os.path.join(options.output_dir,
    "%s.coords" % input_root)
  options.overlaps_path = options.group_coords_path.replace(".coords",
    ".olap")
  options.use_existing_group_coords = False
  if (options.use_existing_overlaps): #{
    if (os.path.isfile(options.group_coords_path)): #{
      options.use_existing_group_coords = True
    #} end if
    if (not os.path.isfile(options.overlaps_path)): #{
      options.use_existing_overlaps = False
    #} end if
  #} end if
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "prediction", options.output_dir)
    if (options.realign): #{
      options.realign_dir = os.path.join(options.output_dir, "realignment")
      CheckDirPath(options.realign_dir, "intermediate realignment",
        path_errors, create=True)
    #} end if
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
    if (CheckPaths(options)): #{
      try: #{
        event_predictor = EventPredictionCls(options)
        WriteCommand(event_predictor, sys.argv)
        event_predictor.PredictEvents()
      except (MyError), e:
        ErrMsg("ERROR while getting biologically typed events in "
               "anomalous contigs:\n  %s" % e)
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
