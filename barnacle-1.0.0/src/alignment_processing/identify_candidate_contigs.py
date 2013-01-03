#! /usr/bin/env python
"""
identify_candidate_contigs.py

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
from utils.log import CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  RunOverlapCode, NormalizeChrID, AddChr, NonStandardChr)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (FileBoxCls, CheckFilePath, CheckDirPath,
  EnsureAbsPath, CleanLine)
from utils.subprocesses import RunCommandFromString, STRING_OUT, STREAM_OUT
from parsers.tokenizer import TokenizerCls
from parsers import psl_parser, exonerate_parser
# use psl_parser.parse(file, filters, noblocks=False noline=False)
# use exonerate_parser.parse(file, filters)
# filters: a dictionary with the following keys relevant
#   bestn:  for each query, return all alignments with rank <= bestn
#   unique: for each query with a single best alignment, return that alignment
#   count:  for each query, return the first count alignments
#   qlen:   for each query, return all alignments with query length >= qlen
#   match:  for each query, return all alignments with >= match% of query
#           bases matching
#   identity: for each query, return all alignments with percent identity >=
#             identity
#   ungapped: for each query, return all alignments with only a single block
from contig_with_alignments import (ContigWithAlignmentsCls,
  ContigAlignsError)
from alignment_functions import (ParseAlignmentFile, FixAlign, CalcOverlap,
  WriteBlockCoords, AlignListString, ShortAlignString, NoAlignmentsError)
from split_candidate import SplitCandidateCls
from gene_overlap import (FeatureOverlapCls, AddOverlapToAlign,
  ChooseBestTranscripts)
from gap_filter import GapFilterError

# CONSTANTS
ES_SUCCESS = 0
ES_SPLIT_FINDER_ERR = 1
ES_OPT_ERR = 2
ES_NO_ALIGNMENTS = 3
ES_PATH_ERR = 4
ES_EXCEPTION = 5
MSG_SUCCESS = "CANDIDATE IDENTIFICATION SUCCESS"
MSG_FAIL = "CANDIDATE IDENTIFICATION FAIL"

# Class to identify interesting split alignments (chimeric RNA, exon
# duplication, local inversions) and gapped alignments (duplications
# and inversions)
class CandidateIdentifierCls: #{
  def __init__(self, options, params): #{
    SetupMainClass(self, options, params=params)
    self.candidate_contigs = []
    self.num_contigs       = 0
    self.num_full_aligns   = 0
    self.num_gapped_aligns = 0
    self.more_than_99 = False
    self.gapped_psl_lines = list()
    self.paths = None
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  # parse the alignment file and examine the contig alignments
  # for interesting contigs
  def CheckContigAlignments(self): #{
    # setup all the paths that will be needed
    self.SetupPaths()
    # filters to use when parsing
    filters = {
      'count':self.params['num_aligns'],
      'identity':self.params['min_identity']
    }
    # parse alignment file
    aligns = ParseAlignmentFile(self.options.input_path, filters,
      self.log_info, ParseFunction=self.parse_aligns)
    # if the get exon info flag is set
    if (self.options.add_gene_annotation): #{
      self.AddGeneAnnotation(aligns)
    #} end if
    # identify candidate contigs
    self.IdentifyCandidateContigs(aligns)
    LogMsg(self, "Done identifying candidate contigs")
    LogMsg(self, "Number of contigs with alignment(s): %i" % self.num_contigs)
    LogMsg(self, "Number of contigs with a single best alignment: %i" %
      self.num_full_aligns)
    LogMsg(self, "Number of alignments: %i" % len(aligns))
    LogMsg(self, "Number of split alignments found: %i" %
      len(self.candidate_contigs))
    if (self.options.check_gap): #{
      msg = "Number of gapped alignments found: "
      if (self.more_than_99): #{
        msg += "at least "
      #} end if
      msg += "%i" % self.num_gapped_aligns
      LogMsg(self, msg)
    #} end if
  #} end def

  def SetupPaths(self): #{
    self.paths = dict()
    input_file_name = os.path.basename(self.options.input_path)
    (root, input_ext) = os.path.splitext(input_file_name)
    self.paths['root'] = root
    # use the input file name extension to determine
    # the tool used for the alignment
    if (".psl" == input_ext): #{
      #align_type = "blat"
      self.parse_aligns = psl_parser.parse
    elif (".exon" == input_ext): #{
      if (self.options.check_gap): #{
        raise CandidateIdentifierError("Cannot use gap filter with "
          "exonerate alignments")
      #} end if
      #align_type = "exon"
      self.parse_aligns = exonerate_parser.parse
    else:
      raise CandidateIdentifierError("Unrecognized input file extension: "
        "\"%s\". " % input_ext + "Should be either \".exon\" or \".psl\"")
    #} end if
    # construct the output paths
    output_root = root
    if (None != self.options.contig_set): #{
      output_root += ".%s" % self.options.contig_set
    #} end if
    #split_out_file_name = "%s.split.candidate.%s" % (output_root, align_type)
    split_out_file_name = "%s.split.candidates" % (output_root)
    self.paths['split_out'] = os.path.join(self.options.output_dir,
      split_out_file_name)
    DebugMsg(self, "Writing split candidates to: %s" % self.paths['split_out'])
    if (self.options.check_gap): #{
      gap_out_file_name = "%s.gap.candidates" % output_root
      self.paths['gap_out'] = os.path.join(self.options.output_dir,
        gap_out_file_name)
      if (not self.options.append and os.path.exists(self.paths['gap_out'])): #{
        os.remove(self.paths['gap_out'])
      #} end if
      DebugMsg(self, "Writing gap candidates to: %s" % self.paths['gap_out'])
      # get the contig sequences path from the input path or the given option
      if (None == self.options.ctg_seq_path): #{
        self.paths['ctg_seq'] = self.options.input_path.replace(
          "output", "input")
        self.paths['ctg_seq'] = self.paths['ctg_seq'].replace(input_ext, ".fa")
      else:
        self.paths['ctg_seq'] = self.options.ctg_seq_path
      #} end if
      # ensure that the contig sequences file exists
      CheckFilePath(self.paths['ctg_seq'], "contig sequence")
    else:
      self.paths['gap_out'] = None
    #} end if
    # construct the counts and psl out paths
    #file_name = "%s.counts.%s" % (output_root, align_type)
    file_name = "%s.counts" % (output_root)
    self.paths['counts'] = os.path.join(self.options.output_dir, file_name)
    #file_name = "%s.%s.psl" % (output_root, align_type)
    file_name = "%s.psl" % (output_root)
    self.paths['psl_out'] = os.path.join(self.options.output_dir, file_name)
  #} end def

  def AddGeneAnnotation(self, aligns): #{
    LogMsg(self, "Getting gene annotation info...")
    # create the gene feature overlap output directory
    overlap_dir = self.CreateOverlapOutputDir()
    # if a gene feature coordinates file was given
    if (None != self.options.gene_coords_path): #{
      DebugMsg(self, "Creating new gene feature overlap file...")
      # setup the alignment coords file (it might already exist)
      self.SetupAlignCoords(aligns, overlap_dir)
      # run the overlap code to get a file with the overlaps
      overlap_stream = self.GetOverlapStream()
    #} end if
    # parse the results of the overlap code
    self.ParseOverlaps(aligns, overlap_stream)
    if (not self.options.keep_coords_file): #{
      os.remove(self.paths['align_coords'])
    #} end def
  #} end def

  def CreateOverlapOutputDir(self): #{
    # remove the input file name
    input_dir = os.path.dirname(self.options.input_path)
    # get the alignment group name
    align_group_dir = os.path.dirname(input_dir)
    align_group = os.path.basename(align_group_dir)
    overlap_dir = os.path.join(self.options.output_dir, "gene_overlap",
      align_group)
    CheckDirPath(overlap_dir, "gene feature overlap", create=True)
    return overlap_dir
  #} end def

  def SetupAlignCoords(self, aligns, overlap_dir): #{
    # determine the alignment coordinates file path
    self.paths['align_coords'] = self.GetIntermediateFilePath(overlap_dir,
      "coords", "CoordsFN")
    # if the alignment coordinates file does not already exist
    if (not os.path.exists(self.paths['align_coords']) or
        0 == os.path.getsize(self.paths['align_coords'])):
      # create the alignment coords file
      self.CreateAlignCoordsFile(aligns)
    #} end if
  #} end def

  def GetIntermediateFilePath(self,
      out_dir, ext, file_descriptor="File path"):
    intermed_file_name = "%s.%s" % (self.paths['root'], ext)
    intermed_file_path = os.path.join(out_dir, intermed_file_name)
    DebugMsg(self, "%s: %s" % (file_descriptor, intermed_file_path))
    return intermed_file_path
  #} end def

  def CreateAlignCoordsFile(self, aligns): #{
    DebugMsg(self, "Creating new alignment coordinates file...")
    # open the alignment coordinates file
    fail_msg = "Cannot open alignment coordinates file"
    align_coords_file = FileBoxCls(self.paths['align_coords'], "w", fail_msg)
    # iterate through the alignments
    for id, align in enumerate(aligns): #{
      align = FixAlign(align)
      # REMINDER: use alignment blocks instead!
      WriteBlockCoords(align, id, align_coords_file, use_chr=True)
      #coord_str = "%s %i %i %i" % (align.target,
      #  min(align.tstart, align.tend), max(align.tstart, align.tend),
      #  id)
      #align_coords_file.write(coord_str + "\n")
    #} end for
    align_coords_file.close()
  #} end def

  def GetOverlapStream(self): #{
    #DebugMsg(self, "Running overlap code...")
    return RunOverlapCode(self.options.gene_coords_path,
      self.paths['align_coords'], STREAM_OUT, get_size=True,
      get_closest=True, dpt=self.options.dpt, log_info=self.log_info)
  #} end def

  def ParseOverlaps(self, aligns, overlap_stream): #{
    LogMsg(self, "Parsing gene feature overlaps...")
    parse_start = time.time()
    # parse the gene feature overlaps result stream
    for overlap_line in overlap_stream: #{
      overlap_line = CleanLine(overlap_line)
      ExtremeDebugMsg(self, "OVERLAP LINE: \"%s\"" % overlap_line)
      tokenizer = TokenizerCls(overlap_line, delimiter=" ",
        log_info=self.log_info)
      try:
        feature_overlap = FeatureOverlapCls(tokenizer)
      except ValueError,e:
        raise CandidateIdentifierError("error parsing overlap line: "
          "%s\n%s" % (overlap_line, e))
      # end try
      if (len(aligns) <= feature_overlap.NumID()): #{
        raise CandidateIdentifierError("invalid alignment id: "
          "%i (max is %i)" % (feature_overlap.NumID(), len(aligns)))
      #} end if
      FixAlign(aligns[feature_overlap.NumID()])
      if ("none" != feature_overlap.overlap_type): #{
        ExtremeDebugMsg(self, "Adding feature %s to align %i" %
          (feature_overlap.id_dict["gene"], feature_overlap.NumID()))
      #} end if
      AddOverlapToAlign(aligns[feature_overlap.NumID()], feature_overlap)
    #} end for
    LogMsg(self, "Time spent parsing gene feature overlaps: %s" %
      TimeSpent(parse_start))
  #} end def

  # identify alignment pairs that look like they could represent
  # interesting split alignments
  def IdentifyCandidateContigs(self, aligns): #{
    # TEMP # ExtremeDebugMsg(self, AlignListString(aligns)

    # open the contig sequences file if using the gap filter
    if (self.options.check_gap): #{
      fail_msg = "Cannot open contig sequence file"
      ctg_seq_file = FileBoxCls(self.paths['ctg_seq'], "r", fail_msg)
    else:
      ctg_seq_file = None
    #} end if

    # iterate over the alignments, grouping them by query (i.e. contig)
    contig_align_index = 0
    while (contig_align_index < len(aligns)): #{
      self.num_contigs += 1
      contig = ContigWithAlignmentsCls(contig_align_index, aligns,
        ctg_seq_file, self.paths['gap_out'], self.options, self.log_info)
      ExtremeDebugMsg(self, "-"*80)
      #DebugMsg(self, "Grouping alignments for "
      #  "%s (contig #%i)..." % (contig.id, self.num_contigs))
      DebugMsg(self, "%i) %s" % (self.num_contigs, contig.id))
      ExtremeDebugMsg(self, "  Contig length: %i" % contig.length)
      #LogMsg(self, "Contig align index: %i" % contig_align_index)

      # Select the alignments to consider for the current contig
      # and check for gapped alignments at the same time
      contig.SelectAlignments(aligns)
      if (contig.single_align_found): #{
        self.num_full_aligns += 1
      #} end if
      if (self.options.check_gap and not contig.perfect_align_found): #{
        contig.CheckGappedAlignments()
        self.gapped_psl_lines.extend(contig.gapped_psl_lines)
        self.num_gapped_aligns += contig.num_gapped_aligns
      #} end if
      if (self.params['check_split'] and
          not self.params['use_quick_chooser']):
        # pare down the alignment groups so that
        # only the best alignments remain
        contig.PareAlignmentGroups()
      #} end if

      #LogMsg(self, "# Gaps Found (Finder): %i" %
      #                         self.num_gapped_aligns)
      contig_align_index += contig.num_aligns_to_contig

      if (0 < len(contig.best_aligns)): #{
        if (self.log_info['debug']): #{
          LogMsg(self, "%i best aligns: %s" %
            (len(contig.best_aligns), contig.id))
          ExtremeDebugMsg(self, AlignListString(contig.best_aligns))
        #} end if
      elif (0 < len(contig.align_groups)): #{
        if (self.log_info['debug']): #{
          ExtremeDebugMsg(self, "-"*40)
          LogMsg(self, "%i align groups: %s" %
            (len(contig.align_groups), contig.id))
          for i, group in enumerate(contig.align_groups): #{
            ExtremeDebugMsg(self, "\n".join(["Group %i" % i,
              "  %i) S:%i E:%i Aligns:%i" % (i, group.ctg_start,
              group.ctg_end, len(group.best_aligns)),
              AlignListString(group.best_aligns)]))
          #} end for
        #} end if
      else: # no best aligns or align groups found
        if (not contig.perfect_align_found and
            not contig.single_align_found): #{
          DebugMsg(self, "No partial aligns selected: %s" % contig.id)
        #} end if
        continue
      #} end if

      # examine pairs of the chosen alignments
      if (self.params['use_quick_chooser']): #{
        self.ExamineBestAlignsPairwise(contig)
      else:
        self.ExamineAlignGroupsPairwise(contig)
      #} end if
    #} end while
    DebugMsg(self, "-"*80)
    # close the contig sequences file if using the gap filter
    if (self.options.check_gap): #{
      ctg_seq_file.close()
    #} end if
  #} end def

  def ExamineBestAlignsPairwise(self, contig): #{
    # iterate over the best alignments
    for a1_index, align1 in enumerate(contig.best_aligns): #{
      #LogMsg(self, "A1: %i" % a1_index) # DEBUG
      if (self.CheckMito(align1)): #{
        continue
      #} end if
      for align2 in contig.best_aligns[a1_index+1:]: #{
        #LogMsg(self, "A1:%i, A2:%i" %
        #  (contig.best_aligns.index(align1),
        #   contig.best_aligns.index(align2))) # DEBUG
        if (self.CheckMito(align2)): #{
          continue
        #} end if
        # examine the alignment pair
        self.ExamineAlignmentPair(align1, align2, contig)
      #} end for
    #} end for
  #} end def

  def ExamineAlignGroupsPairwise(self, contig): #{
    if (2 > len(contig.align_groups)): #{
      ExtremeDebugMsg(self, "No alignment group pairs to compare")
      return
    #} end if
    ExtremeDebugMsg(self, "-"*40 + "\nPairwise comparing %s alignment "
      "groups..." % len(contig.align_groups))
    # iterate over the pairs of alignment groups
    for g1_index, group1 in enumerate(contig.align_groups): #{
      #ExtremeDebugMsg(self, "G1: %i" % g1_index)
      for group2 in contig.align_groups[g1_index+1:]: #{
        #ExtremeDebugMsg(self, "G1:%i, G2:%i" %
        #  (contig.align_groups.index(group1),
        #   contig.align_groups.index(group2)))
        # REMINDER: if the alignments in the groups cannot possibly
        #   represent enough of the contig, skip the whole group pair
        if (self.GroupPairNoGood(group1, group2, contig.length)): #{
          ExtremeDebugMsg(self, "  Pair no good")
          continue
        #} end if
        # iterate over the alignments in the groups
        for align1 in group1.best_aligns: #{
          if (self.CheckMito(align1)): #{
            ExtremeDebugMsg(self, "Skipping mitochondrial alignment")
            continue
          #} end if
          for align2 in group2.best_aligns: #{
            if (self.CheckMito(align2)): #{
              ExtremeDebugMsg(self, "Skipping mitochondrial alignment")
              continue
            #} end if
            # examine the alignment pair
            self.ExamineAlignmentPair(align1, align2, contig)
          #} end for
        #} end for
      #} end for
    #} end for
  #} end def

  def GroupPairNoGood(self, group1, group2, ctg_len): #{
    ExtremeDebugMsg(self, "Checking whether groups could have high "
      "enough contig representation.")
    # REMINDER: implement calculation for best_possible_contig_rep_fraction
    overlap = max(0, CalcOverlap(group1.ctg_start, group1.ctg_end,
      group2.ctg_start, group2.ctg_end)[0])
    best_possible_contig_rep = group1.Span() + group2.Span() - overlap
    best_possible_contig_rep_fraction = (
      float(best_possible_contig_rep) / float(ctg_len))
    if (self.log_info['debug']): #{
      data = list([
        "Group A: %i-%i=%i" %
          (group1.ctg_start, group1.ctg_end, group1.Span()),
        "Group B: %i-%i=%i" %
          (group2.ctg_start, group2.ctg_end, group2.Span()),
        "Overlap: %i" % (overlap),
        "Best Rep: %i" % (best_possible_contig_rep),
        "Ctg Len: %i" % (ctg_len),
        "Best Fraction: %.3f" % (best_possible_contig_rep_fraction)
      ])
      ExtremeDebugMsg(self, ", ".join(data))
    #} end if
    # check whether the contig representation could possibly be high enough
    if (self.params['min_ctg_represented'] >
        best_possible_contig_rep_fraction):
      ExtremeDebugMsg(self, "Group pair could not possibly represent enough "
        "contig, skipping group pair")
      return True
    #} end if
    return False
  #} end def

  # skip alignments to mitochondrial DNA if the no-mito option was used
  def CheckMito(self, align): #{
    if (self.params['discard_mito'] and
        #align.target in ["M","MT","chrM","chrMT"]): #{}
        "M" == NormalizeChrID(align.target)): #{
      ExtremeDebugMsg(self, "Skipping mitochondrial alignment...")
      return True
    #} end if
    return False
  #} end def

  def ExamineAlignmentPair(self, align1, align2, contig): #{
    # create a potential candidate contig from the two alignments
    candidate_contig = SplitCandidateCls(align1, align2,
      self.options.add_gene_annotation, contig.gapped_event_found)
    ExtremeDebugMsg(self, "\n".join(["-"*10,
      "A1: %s" % ShortAlignString(align1), "A2: %s" % ShortAlignString(align2),
      "Contig Represented: %i Contig Represented Fraction: %.3f" %
      (candidate_contig.contig_rep, candidate_contig.contig_rep_fraction)]))
    # if the alignments should have been merged,
    # do not consider the pair
    if (self.params['min_merge_overlap'] <=
        candidate_contig.ctg_overlap_fraction):
      ExtremeDebugMsg(self, "\n".join([
        "Skipping pair that should have been merged",
        ShortAlignString(align1), ShortAlignString(align2),
        "Ctg overlap fraction: %.3f" %
        candidate_contig.ctg_overlap_fraction]))
      return
    #} end if
    # check whether the contig representation is high enough
    if (self.params['min_ctg_represented'] >
        candidate_contig.contig_rep_fraction):
      ExtremeDebugMsg(self, "Not enough contig represented, "
        "skipping alignment pair")
      return
    #} end if
    # if we get to this point, the alignment pair is worth keeping
    ExtremeDebugMsg(self, "Good split candidate identified!")
    # label the event topology
    candidate_contig.LabelTopology(self.options.min_end_dup_fract)
    # check the multi-mapping status of the alignments
    candidate_contig.CheckMultiMapping()
    if (self.options.add_gene_annotation): #{
      # choose the best guesses for transcript(s) the contig represents
      ChooseBestTranscripts(candidate_contig.align1,
        self.options.transcript_selection_buffer)
      ChooseBestTranscripts(candidate_contig.align2,
        self.options.transcript_selection_buffer)
    #} end if
    # skip alignments involving non-standard chromosomes
    if (NonStandardChr(candidate_contig.align1.target) or
        NonStandardChr(candidate_contig.align2.target)): #{
      DebugMsg(self, "Skipping non-standard chromosome for contig %s: %s/%s" %
        (contig.id, candidate_contig.align1.target,
        candidate_contig.align2.target))
      return
    #} end if
    # add the candidate contig to the list
    self.candidate_contigs.append(candidate_contig)
    #LogMsg(self, "-----") # DEBUG
  #} end def

  def WriteCounts(self): #{
    # open the counts file
    fail_msg = "Cannot open counts file"
    counts_file = FileBoxCls(self.paths['counts'], "w", fail_msg)
    # write the number of split alignments found
    counts_file.WriteLine("Split: %i" % len(self.candidate_contigs))
    # if gapped alignments were also checked for
    if (self.options.check_gap): #{
      # write the number of gapped alignments found
      msg = "Gapped: "
      if (self.more_than_99): #{
        msg += "at least "
      #} end if
      msg += "%i" % self.num_gapped_aligns
      counts_file.WriteLine(msg)
    #} end if
    counts_file.WriteLine("COMPLETE")
    # close the counts file
    counts_file.close()
  #} end def

  #Generate the output file
  def Output(self, append): #{
    # open the output file in the appropriate mode
    if append: #{
      mode = "a"
    else:
      mode = "w"
    #} end if
    fail_msg = "Cannot open split alignment output file"
    out = FileBoxCls(self.paths['split_out'], mode, fail_msg)
    if (self.params['output_psl']): #{
      if (self.candidate_contigs[0].align1.method == "blat"): #{
        fail_msg = "Cannot open alignment psl output file"
        psl_out = FileBoxCls(self.paths['psl_out'], mode, fail_msg)
        DebugMsg(self,
          "Writing alignment lines to %s" % self.paths['psl_out'])
        # write out the alignment lines for the gapped alignment events found
        for psl_line in self.gapped_psl_lines: #{
          psl_out.Write(psl_line)
        #} end for
      else:
        # only write out psl lines for blat alignments
        self.params['output_psl'] = False
      #} end if
    #} end if
    # write the split alignment details to the output file
    for candidate_contig in self.candidate_contigs: #{
      # skip non-standard chromosomes
      #chr_patt = r"\A(chr)?(\d+|[XY]|MT?)\Z"
      #if (None == re.search(chr_patt, candidate_contig.align1.target) or
      #    None == re.search(chr_patt, candidate_contig.align2.target)):
      #if (NonStandardChr(candidate_contig.align1.target) or
      #    NonStandardChr(candidate_contig.align2.target)): #{
      #  DebugMsg(self, "Skipping non-standard chromosome: %s/%s" %
      #    (candidate_contig.align1.target, candidate_contig.align2.target))
      #  continue
      #} end if
      #if ("chr" != candidate_contig.align1.target[0:3]): #{
      #  candidate_contig.align1.target = ("chr%s" %
      #    candidate_contig.align1.target)
        #LogMsg(self, "  Target: %s" %
        #                        candidate_contig.align1.target)
        #msg = ("Improperly formatted alignment: %s" %
        #       candidate_contig.Details())
        #raise CandidateIdentifierError(msg)
      #} end if
      #if ("chr" != candidate_contig.align2.target[0:3]): #{
      #  candidate_contig.align2.target = ("chr%s" %
      #    candidate_contig.align2.target)
      #} end if
      candidate_contig.align1.target = AddChr(candidate_contig.align1.target)
      candidate_contig.align2.target = AddChr(candidate_contig.align2.target)
      ExtremeDebugMsg(self, "Writing line to %s:\n  %s" %
        (out.path, candidate_contig.Details()))
      out.WriteLine(candidate_contig.Details())
      if (self.params['output_psl']): #{
        psl_out.Write(candidate_contig.align1.psl())
        psl_out.Write(candidate_contig.align2.psl())
      #} end if
    #} end for
    out.close()
    if (self.params['output_psl']): #{
      psl_out.close()
    #} end if
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class CandidateIdentifierError(MyError): #{
  """Exception raised for errors encountered by the CandidateIdentifierCls class"""
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Examines contig alignments to find contigs "
    "resulting in interesting split or gapped alignments potentially "
    "representing chimeric transcripts.")
  args = [ "ALIGNMENT_FILE", "OUTPUT_DIR", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--split-candidates",
                    action="store_true", dest="check_split",
                    help="Check alignments for contigs matching split "
                         "candidate alignment signatures. [default]")
  parser.add_option("--no-split-candidates",
                    action="store_false", dest="check_split",
                    help="Do not check alignments for contigs matching split "
                         "candidate alignment signatures.")
  parser.add_option("--gap-candidates",
                    action="store_true", dest="check_gap",
                    help="Check alignments for contigs matching gapped "
                         "candidate alignment signatures. [default]")
  parser.add_option("--no-gap-candidates",
                    action="store_false", dest="check_gap",
                    help="Do not check alignments for contigs matching "
                         "gapped candidate alignment signatures.")
  alignment_parsing_grp = OptionGroup(parser, "Alignment Parsing Options")
  alignment_parsing_grp.add_option("-n", "--num-aligns",
                    type='int', dest="num_aligns", metavar="N",
                    help="For each contig extract from the alignments file "
                         "the N highest scoring alignments that pass the PID "
                         "filtering criteria used. [default: %default]")
  alignment_parsing_grp.add_option("-i", "--min-identity",
                    type='float', dest="min_identity", metavar="PID",
                    help="Extract from the alignments file only alignments "
                         "with a percent identity of at least PID. "
                         "[default: %default]")
  alignment_parsing_grp.add_option("--genes",
                    action="store_true", dest="add_gene_annotation",
                    help="Record the gene features overlapped by each "
                         "alignment [default]")
  alignment_parsing_grp.add_option("--no-genes",
                    action="store_false", dest="add_gene_annotation",
                    help="Do not record the gene features overlapped by any "
                         "alignments")
  alignment_parsing_grp.add_option("--gene-coords",
                    metavar="FILE", dest="gene_coords_path",
                    help="Determine whether the alignments overlap any gene "
                         "features using the coordinates in FILE")
  alignment_parsing_grp.add_option("--transcript-selection-buffer",
                    type="int", metavar="N",
                    help="When selecting which transcript(s) a contig "
                         "represents, select all transcripts that overlap "
                         "the alignments by amounts at most Nbp less than "
                         "the largest amount of overlap. [default: %default]")
  alignment_parsing_grp.add_option("--keep-coords-file",
                    action="store_true",
                    help="After gene feature-overlap processing, do not remove "
                         "the alignment-coordinates file that was produced.")
  parser.add_option_group(alignment_parsing_grp)
  alignment_grouping_grp = OptionGroup(parser, "Alignment Grouping Options")
  alignment_grouping_grp.add_option("--single-align",
                    type='float', dest="longest_single_align",
                    metavar="FRACTION",
                    help="If any single alignment for a contig uses at least "
                         "FRACTION of the total contig length, do not look "
                         "for any split or gapped alignments for that contig. "
                         "[default: %default]")
  alignment_grouping_grp.add_option("--merge-overlap",
                    type='float', dest="min_merge_overlap",
                    metavar="FRACTION",
                    help="If two alignments represent portions of the contig "
                         "that overlap by at least FRACTION of the length of "
                         "the shorter alignment, group the alignments "
                         "together. [default: %default]")
  parser.add_option_group(alignment_grouping_grp)
  alignment_selecting_grp = OptionGroup(parser, "Alignment Selection Options")
  alignment_selecting_grp.add_option("--smart-chooser",
                    action="store_false", dest="use_quick_chooser",
                    help="Use the smarter, but slower, choose best "
                         "alignments method [default]")
  alignment_selecting_grp.add_option("--quick-chooser",
                    action="store_true", dest="use_quick_chooser",
                    help="Use the simple, but faster, choose best alignments "
                         "method")
  alignment_selecting_grp.add_option("--maintain-pared-groups",
                    action="store_true", dest="maintain_pared_groups",
                    help="Use paring info to only add an alignment to a "
                         "group if the alignment will not later be pared; "
                         "only relevant when using smart chooser [default]")
  alignment_selecting_grp.add_option("--pare-after-grouping",
                    action="store_false", dest="maintain_pared_groups",
                    help="Do not pare alignments until all alignments for "
                         "the current contig have been grouped; only "
                         "relevant when using smart chooser")
  alignment_selecting_grp.add_option("--mm-min-score",
                    type='float', metavar="F", dest="mm_min_score_fract",
                    help="When deciding whether an alignment group "
                         "multi-maps, only use alignments with scores at "
                         "least F of the maximum score in the group "
                         "[default: %default]")
  alignment_selecting_grp.add_option("--mm-max-pid-diff",
                    type='float', metavar="F", dest="mm_max_pid_diff",
                    help="When deciding whether an alignment group "
                         "multi-maps, only use alignments with percent "
                         "identities at most F less than the maximum percent "
                         "identify in the group; only relevant when using "
                         "smart chooser [default: %default]")
  alignment_selecting_grp.add_option("--mm-max-pid-diff-gap",
                    type='float', metavar="F", dest="mm_max_pid_diff_gap",
                    help="When deciding whether the alignment group used to "
                         "look for gap candidates multi-maps, only use "
                         "alignments with percent identities at most F less "
                         "than the maximum percent identify in the group; "
                         "only relevant when using smart chooser [default: "
                         "%default]")
  alignment_selecting_grp.add_option("--min-score",
                    type='float', metavar="F", dest="min_score_fract",
                    help="When paring alignment groups, only keep alignments "
                         "with scores at least F of the maximum score "
                         "in the group [default: %default]")
  alignment_selecting_grp.add_option("--max-pid-diff",
                    type='float', metavar="F", dest="max_pid_diff",
                    help="When paring alignment groups, only keep alignments "
                         "with percent identities at most F less than the "
                         "maximum percent identify in the group; only "
                         "relevant when using smart chooser [default: "
                         "%default]")
  alignment_selecting_grp.add_option("--prefer-exons",
                    action="store_true",
                    help="When paring alignment groups, if any alignments in "
                         "the group overlap exons, discard all alignments "
                         "in the group that do not overlap any exon")
  alignment_selecting_grp.add_option("--prefer-spliced",
                    action="store_true",
                    help="When paring alignment groups, if any alignments in "
                         "the group are spliced, discard all unspliced "
                         "alignments in the group; only relevant when using "
                         "smart chooser")
  parser.add_option_group(alignment_selecting_grp)
  split_candidate_grp = OptionGroup(parser, "Split Candidate Options")
  split_candidate_grp.add_option("--ctg-rep",
                    type='float', dest="min_ctg_represented",
                    metavar="FRACTION",
                    help="Only report split alignments when the total amount "
                         "of the contig involved in the two alignments are at "
                         "least FRACTION of the total contig length. "
                         "[default: %default]")
  split_candidate_grp.add_option("--min-end-dup-fract",
                    type='float', dest="min_end_dup_fract",
                    metavar="F",
                    help="When labeling alignment topologies, use the "
                         "end-duplication label if the target regions "
                         "overlap by at least F of the smaller target "
                         "region [default: %default]")
  parser.add_option_group(split_candidate_grp)
  gap_realignment_grp = OptionGroup(parser, "Gap Realignment Options")
  gap_realignment_grp.add_option("--ctg-file",
                    dest="ctg_seq_path", metavar="FILE",
                    help="Read in contig sequences from FILE. Otherwise, the "
                         "contig sequences file path is inferred from the "
                         "alignments path. This option is ignored without "
                         "the --use-gap option.")
  gap_realignment_grp.add_option("--gap-realigner",
                    metavar="FILE",
                    help="Use FILE as the program to realign gap sequence of "
                         "gapped alignments")
  gap_realignment_grp.add_option("--gap-config",
                    metavar="FILE",
                    help="FILE is the configuration file for the gap "
                         "realigner.")
  gap_realignment_grp.add_option("--gap-min-size",
                    type='int', dest="min_gap_size", metavar="N",
                    help="When processing alignments, only realign gaps of "
                         "at least N bases [default: %default]")
  gap_realignment_grp.add_option("--gap-check-min-pid",
                    type='float', dest="min_gap_check_pid", metavar="F",
                    help="When processing alignments, only realign gaps "
                         "when the initial alignment has a percent identity "
                         "of at least F. [default: %default]")
  gap_realignment_grp.add_option("--gap-check-min-len",
                    type='float', dest="min_gap_check_len_fract", metavar="F",
                    help="When processing alignments, only realign gaps "
                         "when the initial alignment involves at least F of "
                         "the total contig length [default: %default]")
  gap_realignment_grp.add_option("--gap-check-min-score",
                    type='float', dest="min_gap_check_score_fract",
                    metavar="F",
                    help="When processing alignments, only realign gaps "
                         "when the initial alignment score is at least  F "
                         "of the total contig length [default: %default]")
  gap_realignment_grp.add_option("--gap-max-num-aligns",
                    type="int", metavar="N",
                    help="If a contig has more than N \"good\" alignments, "
                         "do not attempt gap realignment for any of them. "
                         "[default: %default]")
  gap_realignment_grp.add_option("--gap-min-identity",
                   type='float', dest="min_gap_pid", metavar="PID",
                   help="When realigning gaps, only report candidates when "
                        "the gap sequence aligns back to the rest of the "
                        "contig with a percent identity of at least PID "
                        "[default: %default]")
  gap_realignment_grp.add_option("--gap-min-fraction",
                   type='float', dest="min_gap_fract", metavar="F",
                   help="When realigning gaps, only report candidates when "
                        "at least F of the bases in the gap sequence align "
                        "back to the rest of the contig [default: %default]")
  gap_realignment_grp.add_option("--gap-max-len",
                  type="int", metavar="N",
                  help="Only run the gap-realigner on contigs shorter "
                    "than Nbp long. [default: %default]")
  gap_realignment_grp.add_option("--gap-debug",
                  action="store_true",
                  help="Use the debug option when calling the gap "
                       "realignment tool.")
  parser.add_option_group(gap_realignment_grp)
  misc_grp = OptionGroup(parser, "Miscellaneous Options")
  misc_grp.add_option("--no-mito",
                    action="store_true", dest="no_mito",
                    help="Do not report results involving mitochondrial DNA. "
                         "[default]")
  misc_grp.add_option("--allow-mito",
                    action="store_false", dest="no_mito",
                    help="Allow results involving mitochondrial DNA.")
  misc_grp.add_option("--output-psl",
                    action="store_true",
                    help="Output psl alignment lines for events found "
                         "[default]")
  misc_grp.add_option("--no-output-psl",
                    action="store_false", dest="output_psl",
                    help="Do not output psl alignment lines for events found")
  misc_grp.add_option("--contig-set",
                    help="Include the contig set name in the output file "
                         "names (e.g. main, adj, etc.)")
  misc_grp.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  misc_grp.add_option("-a", "--append",
                    action="store_true", dest="append",
                    help="Append to the output file, rather than "
                         "overwriting it.")
  misc_grp.add_option("--log-file",
                    dest="log_file_name", metavar="FILE",
                    help="Log all messages in FILE")
  misc_grp.add_option("-d", "--debug",
                    action="store_true", dest="debug",
                    help="Print debug information while the program runs.")
  misc_grp.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_grp)
  parser.set_defaults(check_split=True,
                      check_gap=True,
                      num_aligns=500,
                      min_identity=40.0,
                      add_gene_annotation=True,
                      transcript_selection_buffer=10,
                      keep_coords_file=False,
                      longest_single_align=0.999,
                      min_merge_overlap=0.8,
                      use_quick_chooser=False,
                      maintain_pared_groups=True,
                      mm_min_score_fract=0.8,
                      mm_max_pid_diff=1.0,
                      mm_max_pid_diff_gap=5.0,
                      min_score_fract=0.85,
                      max_pid_diff=0.5,
                      prefer_exons=False,
                      prefer_spliced=False,
                      min_ctg_represented=0.85,
                      min_end_dup_fract=0.80,
                      min_gap_size=4,
                      min_gap_check_pid=40.0,
                      min_gap_check_len_fract=0.4,
                      min_gap_check_score_fract=0.4,
                      gap_max_num_aligns=3,
                      min_gap_pid=0.95,
                      min_gap_fract=0.3,
                      gap_max_len=50000,
                      gap_debug=False,
                      no_mito=True,
                      output_psl=True,
                      dpt=False,
                      append=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  path_errors = list()
  CheckFilePath(options.input_path, "contig alignment", path_errors)
  CheckDirPath(options.output_dir, "output", path_errors, create=True)
  opts_good = True
  if (options.check_gap): #{
    # check that the config file path was given and exists
    if (None == options.gap_config): #{
      ErrMsg("you must specify the path to the gap realigner config file "
             "when using the gap filter (--gap-config option)")
      opts_good = False
    else:
      CheckFilePath(
        options.gap_config, "gap realigner tool config", path_errors)
    #} end if
    # check that the gap realigner tool path was given and
    # that the tool exists
    if (None == options.gap_realigner): #{
      ErrMsg("you must specify the path to the gap realigner tool when "
             "using the gap filter (--gap-realigner option)")
      opts_good = False
    else:
      CheckFilePath(options.gap_realigner, "gap realigner tool", path_errors)
    #} end if
  #} end if
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # the paths are good if there are no path errors and no conflicting options
  return (opts_good and 0 == len(path_errors))
#} end def

def GetParams(options): #{
  params = {}
  params['longest_single_align']   = options.longest_single_align
  params['min_merge_overlap']      = options.min_merge_overlap
  params['min_ctg_represented']    = options.min_ctg_represented
  params['num_aligns']             = options.num_aligns
  params['min_identity']           = options.min_identity
  params['discard_mito']           = options.no_mito
  params['check_split']            = options.check_split
  params['gap_realign_path']       = options.gap_realigner
  params['min_gap_size']           = options.min_gap_size
  params['min_gap_pid']            = options.min_gap_pid
  params['min_gap_check_pid']      = options.min_gap_check_pid
  # if using the gap filter
  if (options.check_gap): #{
    # get the argument character limit
    output = RunCommandFromString("getconf ARG_MAX", stdout=STRING_OUT,
      dpt=options.dpt)[0]
    params['arg_max'] = int(output)
    options.arg_max = int(output)
    params['gap_config'] = options.gap_config
  #} end if
  params['use_blat_server'] = False
  if (options.prefer_exons): #{
    options.add_gene_annotation = True
  #} end if
  params['use_quick_chooser']     = options.use_quick_chooser
  params['maintain_pared_groups'] = options.maintain_pared_groups
  params['output_psl']            = options.output_psl
  params['debug']                 = options.debug
  return params
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.input_path = EnsureAbsPath(args[0])
    options.output_dir = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      params = GetParams(options)
      try:
        candidate_indentifier = CandidateIdentifierCls(options, params)
        WriteCommand(candidate_indentifier, sys.argv)
        candidate_indentifier.CheckContigAlignments()
      except NoAlignmentsError, e:
        ErrMsg(e)
        return ES_NO_ALIGNMENTS
      except (MyError), e:
        ErrMsg("Error while searching for candidate contigs: %s" % e)
        return ES_SPLIT_FINDER_ERR
      # end try
      if (0 < len(candidate_indentifier.candidate_contigs)): #{
        candidate_indentifier.Output(append=options.append)
      #} end if
      if(options.debug): #{
        ErrMsg("Output to: %s" % options.output_dir)
      #} end if
      candidate_indentifier.WriteCounts()
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify the path to a contig-to-genome alignment "
      "file [ALIGNMENT_FILE] and an output directory path [OUTPUT_DIR]")
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
