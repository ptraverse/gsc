#! /usr/bin/env python
"""
grouped_candidate.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# TODO

# import standard modules
import re

# import custom modules
from utils.error import MyError
from utils.general import (IntOrNAString, IsEmpty, IsUTR, IsNotUTR,
  IsNonCoding, IsNotNonCoding, RemoveExonTypeLabel, LEFT, RIGHT)
from utils.messages import ErrMsg
from utils.multi_dict import MultiDictCls
from parsers.tokenizer import (TokenizerCls, ParseWarningCls, GetFieldAndValue,
  GetFieldValue)
from alignment_processing.alignment_functions import CalcOverlap
from breakpoint import BreakpointCls, BreakpointPairString
from coord_pair import CoordPairCls, CoordPairError, BetweenCoords
from annotation.create_gene_feature_coords import INTRON

# CONSTANTS
START = "start"
END   = "end"

# topologies:


class GroupedCandidateCls: #{
  def __init__(self, data_str, check_data=False): #{
    self.SetupParseFunctions()
    self.tail         = None
    self.score        = None
    self.warnings     = list()
    self.check_data   = False
    #self.meta_match   = None
    self.meta_fields  = dict()
    self.fail_reasons = set()
    try:
      self.ParseDataString(data_str, check_data)
    except (CoordPairError, ValueError), e:
      raise GroupedCandidateError("cannot parse data string: "
        "%s\n%s" % (e, data_str))
    # end try
    self.parts = (START, END)
  #} end def

  def SetupParseFunctions(self): #{
    self.ParseField = dict({
      'CONTIG':                 self.ParseContigInfo,
      'ALIGN_A':                self.ParseAlignmentA,
      'ALIGN_B':                self.ParseAlignmentB,
      'ALIGNER':                self.ParseAlignerUsed,
      'READ_TO_CTG':            self.ParseReadToContigAll,
      'READ_TO_CTG_UNIQUE':     self.ParseReadCoverage,
      'AVG_READ_TO_CTG':        self.ParseAvgReadToContigAll,
      'AVG_READ_TO_CTG_UNIQUE': self.ParseAvgReadCoverage,
      'NUM_GROUPS':             self.ParseNumGroups,
      'OVERLAPPING_GENES':      self.ParseMainGenes,
      'NEARBY_GENES':           self.ParseNearbyGenes,
      'BREAKPOINT_GENES':       self.ParseBPGenes,
      'ENDS':                   self.ParseEnds,
      'REPEATS':                self.ParseRepeats,
      'JUNCTION':               self.ParseBreakpointPair,
      'BREAKPOINTS':            self.ParseBreakpointPair,
      'BLOCKS':                 self.ParseBlocks,
      'EVENT_SEQ':              self.ParseEventSequence,
      'INS_SEQ':                self.ParseInsertedSequence,
      'META':                   self.ParseMetaData,
      'SCORE':                  self.ParseScore,
      'STATUS':                 self.ParseStatus,
    })
  #} end def

  def SetupParseSplitMetaFunctions(self): #{
    self.ParseMetaField = dict({
      'CO': self.ParseContigOverlap,
      'CR': self.ParseContigRepresented,
      'JD': self.ParseBreakpointDistance,
      'BD': self.ParseBreakpointDistance,
      'TD': self.ParseTargetDistance,
      'QG': self.ParseQueryGapCount,
      'SP': self.ParseSplicedCount,
      'GF': self.ParseGapFoundFlag,
      'MM': self.ParseMultiMapping,
      'EB': self.ParseExonBounds,
    })
  #} end def

  def SetupParseGapMetaFunctions(self): #{
    self.ParseMetaField = dict({
      'CS':     self.ParseContigStrand,
      'ES':     self.ParseEventStrand,
      'GAP':    self.ParseGapCoords,
      'NO_GAP': self.ParseNonGapCoords,
      'DUP':    self.ParseDupCoords,
      'DIST':   self.ParseDupDist,
      'TRIM':   self.ParseTrim,
      'FORM':   self.ParseForm,
      'MM':     self.ParseMultiMapping,
      'EB':     self.ParseExonBounds,
    })
  #} end def

  def __lt__(self, other): #{
    return self.score < other.score
  #} end def

  def ParseDataString(self, data_str, check_data=False): #{
    self.check_data = check_data
    # tokenize the data string
    tokens = TokenizerCls(data_str, " ")
    # parse the candidate ID
    ParseFullEventID(self, tokens.next())
    # parse the topology (to set whether this is a gap candidate or not)
    self.ParseTopology(tokens.next())
    # parse the rest of the data string
    extra_data = list()
    for token in tokens: #{
      (field, value) = GetFieldAndValue(token)
      # if the field is a standard field
      if (field in self.ParseField): #{
        self.ParseField[field](value)
      else:
        # it might be the meta-data field
        if (not self.ParseMetaData(token)): #{
          # or it might just be an extra field
          extra_data.append(token)
        #} end if
      #} end if
    #} end for
    if (0 < len(extra_data)): #{
      self.tail = " ".join(extra_data)
    #} end if
  #} end def

  def ParseTopology(self, field_str): #{
    (field, topology) = GetFieldAndValue(field_str)
    if ("TOPOLOGY" != field): #{
      raise GroupedCandidateError(
        "invalid topology field name: \"%s\"" % field)
    #} end if
    self.topology = topology
    if (self.topology.startswith("gap")): #{
      self.gap = True
      self.SetupParseGapMetaFunctions()
    else:
      self.gap = False
      self.SetupParseSplitMetaFunctions()
    #} end if
  #} end def

  def ParseContigInfo(self, contig_str): #{
    self.contig_info = ContigInfoCls(contig_str)
  #} end def

  def ParseAlignmentA(self, alignment_str): #{
    self.align_info_A = AlignmentInfoCls(alignment_str)
    if (self.align_info_A.ctg_start >= self.align_info_A.ctg_end): #{
      self.AddWarning("incorrectly ordered contig coordinates", alignment_str)
    #} end if
    if (not self.gap and hasattr(self, "align_info_B") and
        self.align_info_A.ctg_start > self.align_info_B.ctg_start):
      self.AddWarning("incorrectly ordered contig coordinates", "%s; %s" %
        (self.align_info_A, self.align_info_B))
    #} end if
    if (hasattr(self, "align_info_B")): #{
      self.alignments = (self.align_info_A, self.align_info_B)
    #} end if
  #} end def

  def ParseAlignmentB(self, alignment_str): #{
    self.align_info_B = AlignmentInfoCls(alignment_str)
    if (self.align_info_B.ctg_start >= self.align_info_B.ctg_end): #{
      self.AddWarning("incorrectly ordered contig coordinates", alignment_str)
    #} end if
    if (not self.gap and hasattr(self, "align_info_A") and
        self.align_info_A.ctg_start > self.align_info_B.ctg_start):
      self.AddWarning("incorrectly ordered contig coordinates", "%s; %s" %
        (self.align_info_A, self.align_info_B))
    #} end if
    if (hasattr(self, "align_info_A")): #{
      self.alignments = (self.align_info_A, self.align_info_B)
    #} end if
  #} end def

  def ParseAlignerUsed(self, aligner_used): #{
    self.aligner = aligner_used
  #} end def

  def ParseReadToContigAll(self, r2c_all_str, require=False): #{
    if ("N/A" == r2c_all_str): #{
      self.read_to_contig_all = None
      if (require): #{
        self.AddWarning("read to contig support not calculated", r2c_all_str)
      #} end if
    else:
      try:
        self.read_to_contig_all = int(r2c_all_str)
      except ValueError, e:
        self.AddWarning("invalid read to contig support value", r2c_all_str,
          "must be an integer")
      # end try
    #} end if
  #} end def

  def ParseAvgReadToContigAll(self, r2c_all_str, require=False): #{
    if ("N/A" == r2c_all_str): #{
      self.avg_read_to_contig_all = None
      if (require): #{
        self.AddWarning("average read to contig support not calculated",
          r2c_all_str)
      #} end if
    else:
      try:
        self.avg_read_to_contig_all = float(r2c_all_str)
      except ValueError, e:
        self.AddWarning("invalid average read to contig support value",
          r2c_all_str, "must be a float")
      # end try
    #} end if
  #} end def

  def ParseReadCoverage(self, read_cov_str, require=False): #{
    if ("N/A" == read_cov_str): #{
      self.read_to_ctg_unique = None
      self.rtc_event_size = None
      if (require): #{
        self.AddWarning("unique read to contig support not calculated",
          read_cov_str)
      #} end if
      return
    #} end if
    try:
      self.read_to_ctg_unique = int(read_cov_str)
      self.rtc_event_size = None
    except ValueError, e:
      read_cov_pattern = r"(?P<num_reads>\d+)\((?P<event_size>\d+)bp:.*\)"
      read_cov_match = re.search(read_cov_pattern, read_cov_str)
      if (None == read_cov_match): #{
        self.read_to_ctg_unique = int(read_cov_str.split("(", 1)[0])
        self.rtc_event_size = -1
        self.AddWarning("could not get unique read coverage", read_cov_str)
      else:
        try:
          self.read_to_ctg_unique = int(read_cov_match.group('num_reads'))
          self.rtc_event_size = int(read_cov_match.group('event_size'))
        except ValueError, e:
          self.AddWarning("invalid read to contig support value",
            read_cov_str, "read support and event size must be integers")
        # end try
      #} end if
    # end try
  #} end def

  def ParseAvgReadCoverage(self, read_cov_str, require=False): #{
    if ("N/A" == read_cov_str): #{
      self.avg_read_to_ctg_unique = None
      if (require): #{
        self.AddWarning("average unique read to contig support "
          "not calculated", read_cov_str)
      #} end if
      return
    #} end if
    try:
      self.avg_read_to_ctg_unique = float(read_cov_str)
    except ValueError, e:
      self.AddWarning("invalid average unique read to contig support value",
        read_cov_str, "must be a float")
    # end try
  #} end def

  def ParseNumGroups(self, num_groups): #{
    self.num_groups = int(num_groups)
  #} end def

  def ParseGenes(self, genes_str, depth, list_end=True): #{
    # create one or two dictionaries of gene information
    #ErrMsg("PARSING GENES: %s" % genes_str)
    genes_str = genes_str.replace(" ","_").replace("UTR","utr")
    #ErrMsg("  after replace: %s" % genes_str)
    if (self.gap): #{
      gene_str_A = genes_str
      gene_str_B = ""
    else:
      try:
        (gene_str_A, gene_str_B) = genes_str.split(";")
      except ValueError, e:
        raise GroupedCandidateError("cannot get gene names from: "
          "%s\n%s" % (genes_str, e))
      # end try
    #} end if
    if ("N/A" == gene_str_A): #{
      gene_dict_A = None
    elif ("none" == gene_str_A.lower() or
        ("(" in gene_str_A and ")" in gene_str_A)):
      gene_dict_A = MultiDictCls(gene_str_A, depth, list_end)
    else:
      # for backwards compatibility
      gene_dict_A = gene_str_A.split(",")
    #} end if
    if ("N/A" == gene_str_B): #{
      gene_dict_B = None
    elif ("none" == gene_str_B.lower() or
        ("(" in gene_str_B and ")" in gene_str_B)):
      gene_dict_B = MultiDictCls(gene_str_B, depth, list_end)
    else:
      # for backwards compatibility
      gene_dict_B = gene_str_B.split(",")
    #} end if
    return (gene_dict_A, gene_dict_B)
    #return (gene_str_A.split(","), gene_str_B.split(","))
  #} end def

  def OldParseGenes(self, genes_str): #{
    genes_str = genes_str.replace(" ","_")
    if (self.gap): #{
      gene_names_A = genes_str
      gene_names_B = ""
    else:
      try:
        (gene_names_A, gene_names_B) = genes_str.split(";")
      except ValueError, e:
        raise GroupedCandidateError("cannot get gene names from: "
          "%s\n%s" % (genes_str, e))
      # end try
    #} end if
    return (gene_names_A.split(","), gene_names_B.split(","))
  #} end def

  def ParseMainGenes(self, genes_str): #{
    # main genes have gene_name(parts_list)
    gene_dicts = self.ParseGenes(genes_str, 2, list_end=True)
    (self.genes_A, self.genes_B) = gene_dicts
    self.main_genes = gene_dicts
  #} end def

  def ParseNearbyGenes(self, genes_str): #{
    # nearby genes have gene_name(transcript(distance))
    gene_dicts = self.ParseGenes(genes_str, 3, list_end=False)
    (self.nearby_A, self.nearby_B) = gene_dicts
    self.nearby_genes = gene_dicts
  #} end def

  def ParseBPGenes(self, genes_str): #{
    # breakpoint genes have gene_name(transcript(part(num_list)))
    gene_dicts = self.ParseGenes(genes_str, 4, list_end=True)
    (self.breakpoint_genes_A, self.breakpoint_genes_B) = gene_dicts
    self.break_genes = gene_dicts
  #} end def

  def ParseEnds(self, ends_str): #{
    (end_A, end_B) = ends_str.split(";")
    if ("N/A" in (end_A, end_B)): #{
      if (end_A != end_B): #{
        self.AddWarning("only one end N/A", ends_str)
      #} end if
      self.ends = (end_A, end_B)
      return
    #} end if
    if (end_A not in ["3'", "5'"]): #{
      self.AddWarning("invalid gene end", end_A,
        "must be either \"3'\" or \"5'\"")
    #} end if
    if (end_B not in ["3'", "5'"]): #{
      self.AddWarning("invalid gene end", end_B,
        "must be either \"3'\" or \"5'\"")
    #} end if
    #if (end_A == end_B): #{
    #  self.AddWarning("invalid gene ends", ends_str,
    #    "cannot be the same")
    #} end if
    self.ends = (end_A[0], end_B[0])
  #} end def

  def ParseRepeats(self, repeats_str): #{
    (repeats_A_str, repeats_B_str) = repeats_str.split(";")
    #if ("None" == repeats_A_str or "N/A" == repeats_A_str): #{
    #  self.repeats_A = [repeats_A_str]
    #else:
    self.repeats_A = repeats_A_str.split(",")
    #} end if
    #if ("None" == repeats_B_str or "" == repeats_B_str): #{
    #  self.repeats_B = list()
    #else:
    self.repeats_B = repeats_B_str.split(",")
    #} end if
    self.repeats = (self.repeats_A, self.repeats_B)
  #} end def

  def ParseBreakpointPair(self, breakpoint_strs): #{
    (breakpoint_strA, breakpoint_strB) = breakpoint_strs.split("-")
    self.breakpointA = BreakpointCls(breakpoint_strA)
    self.breakpointB = BreakpointCls(breakpoint_strB)
    self.breakpoints = (self.breakpointA, self.breakpointB)
  #} end def

  def ParseEventSequence(self, event_sequence_str): #{
    self.event_seq = event_sequence_str.upper()
    if (len(event_sequence_str) != self.align_info_B.ContigSpan()): #{
      self.AddWarning("event sequence length does not match event's "
        "contig span", "%s %s (%ibp)" % (self.align_info_B.ToString(),
        event_sequence_str, len(event_sequence_str)))
    #} end if
  #} end def

  def ParseBlocks(self, blocks_str): #{
    #ErrMsg("Parsing blocks string: %s" % blocks_str)
    (blocks_str_A, blocks_str_B) = blocks_str.split(";")
    self.blocks_A = list()
    self.ParseBlockList(self.blocks_A, blocks_str_A)
    self.blocks_B = list()
    self.ParseBlockList(self.blocks_B, blocks_str_B)
    self.blocks = (self.blocks_A, self.blocks_B)
  #} end def

  def ParseBlockList(self, block_list, blocks_str): #{
    block_coords_warnings = list()
    for block_str in blocks_str.split(","): #{
      #ErrMsg("  BLOCK: %s" % block_str)
      (block_start, block_end) = ParseRegionCoords(block_str,
        "block", block_coords_warnings, allow_single=True)
      #ErrMsg("    Start: %i, End: %i" % (block_start, block_end))
      block_list.append((block_start, block_end))
    #} end def
    if (0 < len(block_coords_warnings)): #{
      self.warnings.extend(block_coords_warnings)
    #} end if
  #} end def

  def ParseInsertedSequence(self, inserted_sequence_str): #{
    self.inserted_seq = inserted_sequence_str
  #} end def

  def ParseMetaData(self, meta_data_str): #{
    self.meta_data = meta_data_str
    meta_tokens = TokenizerCls(meta_data_str, ",")
    extra_meta_data = list()
    meta_field_found = False
    for token in meta_tokens: #{
      try:
        (field, value) = GetFieldAndValue(token)
      except ValueError, e:
        raise GroupedCandidateError("cannot parse meta data: "
          "\"%s\", %s" % (token, e))
      # end try
      # check for erroneous MM field
      if (";MM:" in value): #{
        (value, multi_map_token) = value.rsplit(";", 1)
        if ("MM" in self.ParseMetaField): #{
          multi_map_value = GetFieldValue(multi_map_token)
          self.ParseMetaField['MM'](multi_map_value)
          meta_field_found = True
        else:
          extra_meta_data.append(multi_map_token)
        #} end if
      #} end if
      # if the field is a standard meta field
      if (field in self.ParseMetaField): #{
        self.ParseMetaField[field](value)
        meta_field_found = True
      else:
        # it might just be an extra field
        extra_meta_data.append(token)
      #} end if
    #} end for
    # if no meta fields were found, this is not the meta-data field
    if (not meta_field_found): #{
      self.meta_data = ""
      return False
    #} end if
    self.meta_fields['extra'] = extra_meta_data
    if (self.gap and self.check_data): #{
      try:
        self.CheckGapMetaData()
      except KeyError, e:
        raise GroupedCandidateError("invalid meta data field name "
          "encountered in candidate %s: %s" % (self.IDString(),
          str(e).replace("KeyError: ", "")) + "\n%s" % meta_data_str)
      # end try
    #} end if
    # TEMP
    #for meta_field in self.meta_fields.iteritems(): #{
    #  print "%s: %s" % meta_field
    #} end for
    #print self.WarningsString()
    # TEMP
    return True
  #} end def

  def ParseContigOverlap(self, overlap_str): #{
    overlap_pattern = (r"(?P<ctg_overlap>-?\d+)bp" +
      r"\((?P<ctg_overlap_fract>1.0+|0.\d+)\)")
    overlap_match = re.search(overlap_pattern, overlap_str)
    if (None == overlap_match): #{
      self.AddWarning("could not parse contig overlap information",
        overlap_str)
      return
    #} end if
    self.meta_fields['ctg_overlap'] = int(overlap_match.group('ctg_overlap'))
    self.meta_fields['ctg_overlap_fract'] = float(
      overlap_match.group('ctg_overlap_fract'))
  #} end def

  def ParseContigRepresented(self, ctg_rep_str): #{
    try:
      ctg_rep = float(ctg_rep_str)
    except ValueError,e:
      self.AddWarning("invalid contig representation value", ctg_rep_str,
        "must be a float")
      return
    # end try
    if (0.0 > ctg_rep or 1.00001 <= ctg_rep): #{
      self.AddWarning("invalid contig representation value", ctg_rep_str,
        "must be between 0.0 and 1.0")
      return
    # end try
    self.meta_fields['ctg_rep'] = ctg_rep
  #} end def

  def ParseBreakpointDistance(self, breakpoint_distance_str): #{
    if (breakpoint_distance_str in ["NA", "N/A"]): #{
      self.meta_fields['breakpoint_dist'] = "N/A"
    else:
      try:
        breakpoint_dist = int(breakpoint_distance_str)
        self.meta_fields['breakpoint_dist'] = breakpoint_dist
      except ValueError, e:
        self.AddWarning("invalid breakpoint distance", breakpoint_distance_str,
          "must be an integer")
      # end try
    #} end if
  #} end def

  def ParseTargetDistance(self, target_distance_str): #{
    if (target_distance_str in ["NA", "N/A"]): #{
      self.meta_fields['target_dist'] = "N/A"
    else:
      try:
        target_dist = int(target_distance_str)
        self.meta_fields['target_dist'] = target_dist
      except ValueError, e:
        self.AddWarning("invalid target distance", target_distance_str,
          "must be an integer")
      # end try
    #} end if
  #} end def

  def ParseQueryGapCount(self, query_gap_str): #{
    if (self.CheckCount("query gap", query_gap_str)): #{
      self.meta_fields['query_gap'] = int(query_gap_str)
    #} end if
  #} end def

  def ParseSplicedCount(self, spliced_str): #{
    if (self.CheckCount("spliced", spliced_str)): #{
      self.meta_fields['spliced'] = int(spliced_str)
    #} end if
  #} end def

  def ParseGapFoundFlag(self, gap_found_str): #{
    if (gap_found_str not in ["Y", "N"]): #{
      self.AddWarning("invalid gap found flag", gap_found_str,
        "must be \"Y\" or \"N\"")
      return
    #} end if
    self.meta_fields['gap_found'] = gap_found_str
  #} end def

  def ParseMultiMapping(self, multi_map_str): #{
    if (self.CheckCount("multi-mapping", multi_map_str)): #{
      self.meta_fields['multi_mapped'] = int(multi_map_str)
    #} end if
  #} end def

  def ParseExonBounds(self, exon_bounds_str): #{
    if (self.CheckCount("exon boundaries", exon_bounds_str)): #{
      self.meta_fields['exon_bounds'] = int(exon_bounds_str)
    #} end if
  #} end def

  def CheckCount(self, count_type, count): #{
    if (count not in ["0", "1", "2"]): #{
      self.AddWarning("invalid %s count" % count_type, count,
        "must be 0, 1, or 2")
      return False
    #} end if
    return True
  #} end def

  def ParseContigStrand(self, strand): #{
    if (self.CheckStrand("contig", strand)): #{
      self.meta_fields['ctg_strand'] = strand
    #} end if
  #} end def

  def ParseEventStrand(self, strand): #{
    if (self.CheckStrand("event", strand)): #{
      self.meta_fields['event_strand'] = strand
    #} end if
  #} end def

  def CheckStrand(self, strand_type, strand): #{
    if ("+" != strand and "-" != strand): #{
      self.AddWarning("invalid %s strand" % strand_type, strand,
        "must be \"+\" or \"-\"")
      return False
    #} end if
    return True
  #} end def

  def ParseGapCoords(self, gap_coords_str): #{
    gap_coords_warnings = list()
    (gap_start, gap_end) = ParseRegionCoords(gap_coords_str,
      "gap", gap_coords_warnings)
    if (0 < len(gap_coords_warnings)): #{
      self.warnings.extend(gap_coords_warnings)
      return
    #} end if
    self.meta_fields['gap_start'] = gap_start
    self.meta_fields['gap_end']   = gap_end
    self.meta_fields['gap_coords']   = CoordPairCls(gap_start, gap_end)
  #} end def

  def ParseNonGapCoords(self, no_gap_coords_str): #{
    no_gap_coords_warnings = list()
    (bg_coords_str, ag_coords_str) = no_gap_coords_str.split(";")
    (bg_start, bg_end) = ParseRegionCoords(bg_coords_str, "before gap",
      no_gap_coords_warnings, ordered=False, allow_na=True)
    (ag_start, ag_end) = ParseRegionCoords(ag_coords_str, "after gap",
      no_gap_coords_warnings, ordered=False, allow_na=True)
    if (0 < len(no_gap_coords_warnings)): #{
      self.warnings.extend(no_gap_coords_warnings)
      return
    #} end if
    self.meta_fields['bg_start'] = bg_start
    self.meta_fields['bg_end']   = bg_end
    self.meta_fields['ag_start'] = ag_start
    self.meta_fields['ag_end']   = ag_end
  #} end def

  def ParseDupCoords(self, dup_coords_str): #{
    dup_coords_warnings = list()
    (dup_start, dup_end) = ParseRegionCoords(dup_coords_str, "duplicate",
      dup_coords_warnings, allow_na=True)
    if (0 < len(dup_coords_warnings)): #{
      self.warnings.extend(dup_coords_warnings)
      return
    #} end if
    self.meta_fields['dup_start'] = dup_start
    self.meta_fields['dup_end']   = dup_end
    if ("duplication" in self.topology): #{
      self.meta_fields['dup_coords']   = CoordPairCls(dup_start, dup_end)
    #} end if
  #} end def

  def ParseDupDist(self, dup_dist_str): #{
    if ("N/A" == dup_dist_str): #{
      self.meta_fields['dist'] = dup_dist_str
      return
    #} end if
    try:
      dup_dist = int(dup_dist_str)
    except ValueError, e:
      self.AddWarning("invalid duplication distance", dup_dist_str,
        "must be integer or \"N/A\"")
      return
    # end try
    self.meta_fields['dist'] = dup_dist
  #} end def

  def ParseTrim(self, trim_str): #{
    self.meta_fields['trim'] = trim_str
  #} end def

  def ParseForm(self, form_str): #{
    self.meta_fields['form'] = form_str
  #} end def

  def CheckGapMetaData(self): #{
    if ('ctg_strand' in self.meta_fields and
        self.meta_fields['ctg_strand'] != self.align_info_A.Strand()):
      self.AddWarning("contig strand does not match order of "
        "genomic coordinates", "CS:%s %s" % (self.meta_fields['ctg_strand'],
        self.align_info_A.ToString()))
    #} end if
    if ('event_strand' in self.meta_fields and
        self.meta_fields['event_strand'] != self.align_info_B.Strand()):
      self.AddWarning("event strand does not match order of "
        "genomic coordinates", "ES:%s %s" % (self.meta_fields['event_strand'],
        self.align_info_B.ToString()))
    #} end if
    if ('gap_start' in self.meta_fields and
        'gap_end'   in self.meta_fields and
        self.meta_fields['gap_start'] >= self.meta_fields['gap_end']):
      self.AddWarning("incorrectly ordered gap coordinates",
        "GAP:%s" % self.GapCoordsString())
    #} end if
    if ('bg_start' in self.meta_fields and
        'gap_end'  in self.meta_fields and
        "N/A" == self.meta_fields['bg_start']):
      if (self.meta_fields['gap_end'] >= self.align_info_A.ctg_start): #{
        self.AddWarning("missing \"before gap\" coordinates",
          "CTG:%s GAP:%s" % (self.align_info_A.ContigCoordsString(),
          self.GapCoordsString()))
      #} end if
    #} end if
    if ('ag_start'  in self.meta_fields and
        'gap_start' in self.meta_fields and
        "N/A" == self.meta_fields['ag_start']):
      if (self.meta_fields['gap_start'] <= self.align_info_A.ctg_end): #{
        self.AddWarning("missing \"after gap\" coordinates",
          "CTG:%s GAP:%s" % (self.align_info_A.ContigCoordsString(),
          self.GapCoordsString()))
      #} end if
    #} end if
    if ('dup_start' in self.meta_fields and
        'dup_end'   in self.meta_fields and
        "N/A" != self.meta_fields['dup_start'] and
        self.meta_fields['dup_start'] >= self.meta_fields['dup_end']):
      self.AddWarning("incorrectly ordered duplicate coordinates",
        "DUP:%s" % self.DupCoordsString())
    #} end if
    if ('dist'      in self.meta_fields and
        'dup_start' in self.meta_fields and
        'dup_end'   in self.meta_fields and
        "N/A" != self.meta_fields['dist']):
      if (self.align_info_B.ctg_start <= self.meta_fields['dup_start']): #{
        dup_dist = (self.meta_fields['dup_start'] -
                    self.align_info_B.ctg_end) - 1
      else:
        dup_dist = (self.align_info_B.ctg_start -
                    self.meta_fields['dup_end']) - 1
      #} end if
      if (dup_dist != self.meta_fields['dist']): #{
        self.AddWarning("incorrect distance between duplicates",
          "EVENT:%s DUP:%s DIST:%s" % (self.align_info_B.ContigCoordsString(),
           self.DupCoordsString(), self.meta_fields['dist']))
      #} end if
    #} end if
  #} end def

  def ParseScore(self, score): #{
    self.score = float(score)
  #} end def

  def ParseStatus(self, status_str): #{
    if (status_str.startswith("FAIL")): #{
      fail_pattern = r"FAIL\((?P<fail_reasons>.+)\)"
      fail_match = re.search(fail_pattern, status_str)
      if (None == fail_match): #{
        self.AddWarning("cannot parse fail reasons", status_str)
        return
      #} end if
      self.fail_reasons = set(fail_match.group('fail_reasons').split(","))
    #} end if
  #} end def

  def DataList(self): #{
    data_fields = [
      "%s)" % self.IDString(),
      "TOPOLOGY:%s" % self.topology,
      "CONTIG:%s" % self.contig_info.ToString(),
      "ALIGN_A:%s" % self.align_info_A.ToString(),
      "ALIGN_B:%s" % self.align_info_B.ToString(),
      "ALIGNER:%s" % self.aligner,
      "READ_TO_CTG:%s" % self.ReadToContigAllString(),
      "READ_TO_CTG_UNIQUE:%s" % self.ReadCoverageString(),
      "NUM_GROUPS:%i" % self.num_groups,
      "OVERLAPPING_GENES:%s" % self.MainGenesString(),
      "NEARBY_GENES:%s" % self.NearbyGenesString(),
      "BREAKPOINT_GENES:%s" % self.BPGenesString(),
    ]
    if (hasattr(self, "ends")): #{
      data_fields.append("ENDS:%s" % self.EndsString())
    #} end if
    data_fields.extend([
      "REPEATS:%s" % self.RepeatsString(),
      "META:%s" % self.MetaDataString(),
      "BREAKPOINTS:%s" % self.BreakpointsString(),
    ])
    if (hasattr(self, "avg_read_to_contig_all") and
        None != self.avg_read_to_contig_all): #{
      data_fields.append("AVG_READ_TO_CTG:%s" %
        self.AvgReadToContigAllString())
    #} end if
    if (hasattr(self, "avg_read_to_ctg_unique") and
        None != self.avg_read_to_ctg_unique):
      data_fields.append("AVG_READ_TO_CTG_UNIQUE:%s" %
        self.AvgReadCoverageString())
    #} end if
    if (hasattr(self, "blocks_A") and hasattr(self, "blocks_B")): #{
      data_fields.append("BLOCKS:%s" % self.BlocksString())
    #} end if
    if (self.gap): #{
      data_fields.append("EVENT_SEQ:%s" % self.event_seq)
      if (hasattr(self, "inserted_seq") and not None == self.inserted_seq): #{
        data_fields.append("INS_SEQ:%s" % self.inserted_seq)
      #} end if
    #} end if
    if (not self.PassedFilters()): #{
      data_fields.append("STATUS:%s" % self.Status())
    #} end if
    if (None != self.tail): #{
      data_fields.append(self.tail)
    #} end if
    if (None != self.score): #{
      data_fields.append("SCORE:%.2f" % self.score)
    #} end if
    return data_fields
  #} end def

  def DataString(self): #{
    data_fields = self.DataList()
    return " ".join(data_fields)
  #} end def

  def ReadableString(self): #{
    data_fields = ["-"*80]
    data_fields.extend(self.DataList()[1:])
    # swap order of CONTIG and TOPOLOGY fields
    ctg = data_fields[2]
    data_fields[2] = data_fields[1]
    data_fields[1] = ctg
    if (None != self.score): #{
      data_fields.insert(2, data_fields[-1])
      del data_fields[-1]
    #} end if
    return "\n".join(data_fields)
  #} end def

  def IDString(self): #{
    return "%i%s" % (self.group_id, self.candidate_id)
  #} end def

  def ReadToContigAllString(self): #{
    if (None == self.read_to_contig_all): #{
      return "N/A"
    else:
      return "%i" % self.read_to_contig_all
    #} end if
  #} end def

  def AvgReadToContigAllString(self): #{
    if (None == self.avg_read_to_contig_all): #{
      return "N/A"
    else:
      return "%.3f" % self.avg_read_to_contig_all
    #} end if
  #} end def

  def ReadCoverageString(self): #{
    if (None == self.read_to_ctg_unique): #{
      return "N/A"
    elif (None == self.rtc_event_size): #{
      return "%i" % self.read_to_ctg_unique
    else:
      return "%i(%ibp:%.2f)" % (self.read_to_ctg_unique,
        self.rtc_event_size, self.NormalizedReadToContig())
    #} end if
  #} end def

  def AvgReadCoverageString(self): #{
    if (None == self.avg_read_to_ctg_unique): #{
      return "N/A"
    else:
      return "%.3f" % self.avg_read_to_ctg_unique
    #} end if
  #} end def

  def NormalizedReadToContig(self): #{
    if (0 < self.rtc_event_size): #{
      return float(self.read_to_ctg_unique) / float(self.rtc_event_size)
    else:
      self.AddWarning("non-positive event region size encountered",
        self.rtc_event_size)
      return 0.0
    #} end if
  #} end def

  def GeneSetString(self, gene_set): #{
    if (None == gene_set): #{
      return "N/A"
    #elif (0 == len(gene_set)): #{
    #elif ((hasattr(gene_set, "top_dict") and 0 == len(gene_set.top_dict)) or
    #    0 == len(gene_set)): #{
    elif (0 == len(gene_set)):
      #return "no_exons"
      return "none"
    elif (isinstance(gene_set, MultiDictCls)):
      #return ",".join(gene_set)
      return gene_set.OutputString(sort=True,sort_key=str.lower)
    else:
      # for backwards compatibility
      return ",".join(gene_set)
    #} end if
  #} end def

  def BothGeneSetsString(self, genes_A, genes_B): #{
    gene_names_A = self.GeneSetString(genes_A)
    if (self.gap): #{
      return gene_names_A
    #} end if
    gene_names_B = self.GeneSetString(genes_B)
    return ";".join([gene_names_A, gene_names_B])
  #} end def

  def MainGenesString(self): #{
    return self.BothGeneSetsString(self.genes_A, self.genes_B)
  #} end def

  def NearbyGenesString(self): #{
    if (not hasattr(self, 'nearby_A')): #{
      self.nearby_A = None
    #} end if
    if (not hasattr(self, 'nearby_B')): #{
      self.nearby_B = None
    #} end if
    return self.BothGeneSetsString(self.nearby_A, self.nearby_B)
  #} end def

  def BPGenesString(self): #{
    if (not hasattr(self, 'breakpoint_genes_A')): #{
      self.breakpoint_genes_A = None
    #} end if
    if (not hasattr(self, 'breakpoint_genes_B')): #{
      self.breakpoint_genes_B = None
    #} end if
    return self.BothGeneSetsString(self.breakpoint_genes_A,
      self.breakpoint_genes_B)
  #} end def

  def EndsString(self): #{
    return "%s';%s'" % self.ends
  #} end def

  def RepeatsString(self): #{
    repeats_A_str = ",".join(self.repeats_A)
    if ("" == repeats_A_str): #{
      repeats_A_str = "None"
    #} end if
    repeats_B_str = ",".join(self.repeats_B)
    if ("" == repeats_B_str): #{
      repeats_B_str = "None"
    #} end if
    return ";".join([repeats_A_str, repeats_B_str])
  #} end def

  def MetaDataString(self): #{
    try:
      if (self.gap): #{
        meta_data_list = [
          "CS:%s"  % self.meta_fields['ctg_strand'],
          "ES:%s"  % self.meta_fields['event_strand'],
          "GAP:%s" % self.GapCoordsString(),
          "NO_GAP:%s;%s" %
            (self.BeforeGapCoordsString(), self.AfterGapCoordsString()),
          "DUP:%s"  % self.DupCoordsString(),
          "DIST:%s" % IntOrNAString(self.meta_fields['dist']),
          "TRIM:%s" % self.meta_fields['trim'],
          "FORM:%s" % self.meta_fields['form'],
        ]
      else:
        meta_data_list = [
          "CO:%s"   % self.ContigOverlapString(),
          "CR:%.3f" % self.meta_fields['ctg_rep'],
          "BD:%s"   % IntOrNAString(self.meta_fields['breakpoint_dist']),
          "TD:%s"   % IntOrNAString(self.meta_fields['target_dist']),
          "QG:%i"   % self.meta_fields['query_gap'],
          "SP:%i"   % self.meta_fields['spliced'],
          "GF:%s"   % self.meta_fields['gap_found'],
        ]
      #} end if
      if ("extra" in self.meta_fields): #{
        meta_data_list.extend(self.meta_fields['extra'])
      #} end if
      if ("multi_mapped" in self.meta_fields): #{
        meta_data_list.append("MM:%s" % self.meta_fields['multi_mapped'])
      #} end if
      if ("exon_bounds" in self.meta_fields): #{
        meta_data_list.append("EB:%s" % self.meta_fields['exon_bounds'])
      #} end if
    except KeyError, e:
      raise GroupedCandidateError("invalid meta data field name encountered "
        "in candidate %s: %s" % (self.IDString(),
        str(e).replace("KeyError: ", "")))
    # end try
    return ",".join(meta_data_list)
  #} end def

  def ContigOverlapString(self): #{
    return "%ibp(%.3f)" % (self.meta_fields['ctg_overlap'],
      self.meta_fields['ctg_overlap_fract'])
  #} end def

  def GapSpan(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("No gap span for split candidates.")
    #} end if
    return (self.meta_fields['gap_end'] -
            self.meta_fields['gap_start']) + 1
  #} end def

  def GapCoordsString(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("No gap coords string for split candidates.")
    #} end if
    return "%s-%s" % (self.meta_fields['gap_start'],
      self.meta_fields['gap_end'])
  #} end def

  # realign_query  = portion of gap
  # realign_target = place gap sequence aligns to
  def NewRealignmentCoords(self, trim=False): #{
    #ErrMsg("    Getting realignment coordinates...")
    #realign_qstart = self.align_info_B.ctg_start
    #realign_qend   = self.align_info_B.ctg_end
    realign_query = CoordPairCls(self.align_info_B.ctg_start,
      self.align_info_B.ctg_end)
    #realign_tstart = realign_tend = None
    realign_target = CoordPairCls()
    if ('dup_start' in self.meta_fields and
        "inversion" not in self.topology): #{
      #realign_tstart = self.meta_fields['dup_start']
      #realign_tend   = self.meta_fields['dup_end']
      realign_target = CoordPairCls(self.meta_fields['dup_start'],
        self.meta_fields['dup_end'])
      #gap_coords = CoordPairCls(self.meta_fields['gap_start'],
      #  self.meta_fields['gap_end'])
      #(olapA, fractionA) = CalcOverlap(self.meta_fields['gap_start'],
      #  self.meta_fields['gap_end'], realign_qstart, realign_qend)
      olapA = realign_query.OverlapAmount(self.meta_fields['gap_coords'])
      #(olapB, fractionB) = CalcOverlap(self.meta_fields['gap_start'],
      #  self.meta_fields['gap_end'], self.meta_fields['dup_start'],
      #  self.meta_fields['dup_end'])
      olapB = realign_target.OverlapAmount(self.meta_fields['gap_coords'])
      ctg_str = "ctg:%i-%i(%i)" % (self.align_info_B.ctg_start,
        self.align_info_B.ctg_end, olapA)
      dup_str = "dup:%i-%i(%i)" % (self.meta_fields['dup_start'],
        self.meta_fields['dup_end'], olapB)
      #ErrMsg("Checking which coordinates to use for edge-gap\n"
      #  "gap:%i-%i %s %s" % (self.meta_fields['gap_start'],
      #   self.meta_fields['gap_end'], ctg_str, dup_str))
      if (olapB > olapA): #{
        #realign_qstart = self.meta_fields['dup_start']
        #realign_qend   = self.meta_fields['dup_end']
        #realign_tstart = self.align_info_B.ctg_start
        #realign_tend   = self.align_info_B.ctg_end
        # swap coordinate pairs
        temp = realign_query
        realign_query = realign_target
        realign_target = temp
        #ErrMsg("Using DUP coordinates for query")
      #else:
      #  ErrMsg("Using CTG coordinates for query")
      #} end if
      # trim off bits of realignment extending beyond original gap
      if (trim): #{
        if (self.meta_fields['ctg_strand'] !=
            self.meta_fields['event_strand']): #{
          raise GroupedCandidateError("Realignment trimming not implemented "
            "for inversion gap events.")
        #} end if
        if (realign_query.min < self.meta_fields['gap_coords'].min): #{
          adjustment = self.meta_fields['gap_coords'].min - realign_query.min
          realign_query.min += adjustment
          realign_target.min += adjustment
        #} end if
        if (realign_query.max > self.meta_fields['gap_coords'].max): #{
          adjustment = realign_query.max - self.meta_fields['gap_coords'].max
          realign_query.max -= adjustment
          realign_target.max -= adjustment
        #} end if
      #} end if
    #} end if
    return (realign_query, realign_target)
  #} end def

  def RealignmentCoords(self, trim=False): #{
    (realign_query, realign_target) = self.NewRealignmentCoords(trim)
    return (realign_query.min, realign_query.max,
      realign_target.min, realign_target.max)
  #} end def

  def DupCoordsString(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("No dup coords string for split candidates.")
    #} end if
    if ("N/A" == self.meta_fields['dup_start'] and
        "N/A" == self.meta_fields['dup_end']):
      return "N/A-N/A"
    #} end if
    return "%i-%i" % (self.meta_fields['dup_start'],
      self.meta_fields['dup_end'])
  #} end def

  def BeforeGapCoordsString(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("No before gap coords string for split "
        "candidates.")
    #} end if
    return "%s-%s" % (self.meta_fields['bg_start'], self.meta_fields['bg_end'])
  #} end def

  def AfterGapCoordsString(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("No after gap coords string for split "
        "candidates.")
    #} end if
    return "%s-%s" % (self.meta_fields['ag_start'], self.meta_fields['ag_end'])
  #} end def

  def BreakpointsString(self, sort=False): #{
    pair = (self.breakpointA, self.breakpointB)
    return BreakpointPairString(pair, sort)
  #} end def

  def BlocksString(self): #{
    blocks_str_A = ",".join(["%i-%i" % block for block in self.blocks_A])
    blocks_str_B = ",".join(["%i-%i" % block for block in self.blocks_B])
    return ";".join([blocks_str_A, blocks_str_B])
  #} end def

  def Status(self): #{
    if (0 < len(self.fail_reasons)): #{
      return "FAIL(%s)" % ",".join(sorted(self.fail_reasons))
    else:
      return "PASS"
    #} end if
  #} end def

  def PassedFilters(self): #{
    if (0 < len(self.fail_reasons)): #{
      return False
    else:
      return True
    #} end if
  #} end def

  def AddWarning(self, warning, data, rule=""): #{
    parse_warning = ParseWarningCls(warning, data, rule)
    self.warnings.append(parse_warning)
  #} end def

  def WarningsString(self, indent=""): #{
    warnings_list = list()
    for parse_warning in self.warnings: #{
      message = ("Warning: %s for candidate %s with contig %s: \"%s\"" %
        (parse_warning.warning, self.IDString(),
         self.contig_info.id, parse_warning.data))
      if ("" != parse_warning.rule): #{
        message += ", %s" % parse_warning.rule
      #} end if
      warnings_list.append(message)
    #} end for
    delim = "%s\n" % indent
    return delim.join(warnings_list)
  #} end def

  def TotalContigRepresented(self): #{
    if (self.gap): #{
      # check that the meta-data could be parsed
      #if (None == self.meta_match): #{
      if (0 == len(self.meta_fields)): #{
        raise GroupedCandidateError("could not parse gap meta data for "
          "candidate %s: %s" % (self.IDString(), self.meta_data))
      #} end if
      # before_gap_span + after_gap_span + aligned_gap_span
      ctg_represented = ((self.BeforeGapSpan() +
                          self.AfterGapSpan()) +
                          self.align_info_B.ContigSpan())
    else:
      # align_A_span + align_B_span - overlap
      ctg_represented = ((self.align_info_A.ContigSpan() +
                          self.align_info_B.ContigSpan()) -
                          self.ContigOverlap())
      #ctg_rep_pattern = r",CR:(?P<ctg_rep>[0-9]+\.[0-9]{2}),"
      #ctg_rep_match = re.search(ctg_rep_pattern, self.meta_data)
      #return float(ctg_rep_match.group('ctg_rep'))
    #} end if
    ctg_represented_fraction = (float(ctg_represented) /
      float(self.contig_info.length))
    return ctg_represented_fraction
  #} end def

  def BeforeGapSpan(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("Cannot calculate before gap span for "
        "split candidates.")
    #} end if
    start = self.align_info_A.ctg_start
    end = self.meta_fields['gap_start'] - 1
    # if the gap is at the very beginning of the contig
    if (end <= start): #{
      return 0
    #} end if
    span = (end - start) + 1
    return span
  #} end def

  def AfterGapSpan(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("Cannot calculate before gap span for "
        "split candidates.")
    #} end if
    start = self.meta_fields['gap_end'] + 1
    end = self.align_info_A.ctg_end
    # if the gap is at the very end of the contig
    if (end <= start): #{
      return 0
    #} end if
    span = (end - start) + 1
    return span
  #} end def

  def ContigOverlap(self): #{
    if (self.gap or
        self.align_info_A.ctg_end < self.align_info_B.ctg_start):
      return 0
    #} end if
    overlap = (self.align_info_A.ctg_end - self.align_info_B.ctg_start) + 1
    return overlap
  #} end def

  def MultiMapped(self): #{
    if ('multi_mapped' in self.meta_fields): #{
      return self.meta_fields['multi_mapped']
    else:
      return 0
    #} end if
  #} end def

  def PolyASuspect(self): #{
    polyA_patt = r"^([ACGT])\1*$"
    if (self.gap and
        not self.GapIsInternal() and
        None != re.search(polyA_patt, self.event_seq)):
      return True
    else:
      return False
    #} end if
  #} end def

  def StartOfContig(self): #{
    return 1
  #} end def

  def EndOfContig(self): #{
    return self.contig_info.length
  #} end def

  def GapIsInternal(self): #{
    if (not self.gap): #{
      raise GroupedCandidateError("GapIsInternal() undefined for non-gap "
        "candidates")
    #} end if
    if (self.StartOfContig() == self.meta_fields['gap_start'] or
        self.EndOfContig()   == self.meta_fields['gap_end']):
      return False
    else:
      return True
    #} end if
  #} end def

  def ContigEventRegion(self): #{
    # if this candidate is a gap
    if (self.gap): #{
      #ErrMsg("    Getting gap event region...")
      # figure out which coordinates to use
      (realign_query, realign_target) = self.NewRealignmentCoords()
      #gap_coords = CoordPairCls(self.meta_fields['gap_start'],
      #  self.meta_fields['gap_end'])
      # if the gap is internal to the contig
      if (self.GapIsInternal()): #{
        #ErrMsg("    Gap is internal")
        if (realign_query.Overlaps(realign_target)): #{
          raise GroupedCandidateError("duplicated sequences should not "
            "overlap: %s" % self.IDString())
        #} end if
        # if this is an inversion
        if ("inversion" in self.topology): #{
          # use the inverted region
          region = realign_query
        # if this is a duplication
        else:
          #ErrMsg("    Using region between duplicated sequences...")
          # use the region between the duplicated sequences
          region = BetweenCoords(realign_query, realign_target)
          # if the duplications are adjacent
          if (0 == region.Span()): #{
            region.ResortCoords()
          #} end if
        #} end if
      # if the gap is at an edge of the contig
      else:
        #ErrMsg("    Gap is at edge")
        # if the gap is at the start of the contig
        if (self.StartOfContig() == self.meta_fields['gap_coords'].min): #{
          #ErrMsg("      Gap is at start. Realign Query: %s, Gap: %s" %
          #  (realign_query, gap_coords))
          #left  = min(realign_end, self.meta_fields['gap_end'])
          #right = max(realign_end, self.meta_fields['gap_end'])
          region = CoordPairCls(realign_query.max, self.meta_fields['gap_coords'].max)
        # if the gap is at the end of the contig
        else:
          #ErrMsg("      Gap is at end. Realign Query: %s, Gap: %s" %
          #  (realign_query, gap_coords))
          #left  = min(realign_start, self.meta_fields['gap_start'])
          #right = max(realign_start, self.meta_fields['gap_start'])
          region = CoordPairCls(self.meta_fields['gap_coords'].min, realign_query.min)
        #} end if
      #} end if
    else:
      #left  = min(self.align_info_A.ctg_end, self.align_info_B.ctg_start)
      #right = max(self.align_info_A.ctg_end, self.align_info_B.ctg_start)
      region = CoordPairCls(self.align_info_A.ctg_end,
        self.align_info_B.ctg_start)
    #} end if
    return region
  #} end def

  def ClearBPGenes(self): #{
    #self.breakpoint_genes_A = list()
    #self.breakpoint_genes_B = list()
    #if (None != self.breakpoint_genes_A): #{
    #  self.breakpoint_genes_A.Initialize()
    #} end if
    #if (None != self.breakpoint_genes_B): #{
    #  self.breakpoint_genes_B.Initialize()
    #} end if
    self.ClearGeneSet("breakpoint_genes_A")
    self.ClearGeneSet("breakpoint_genes_B")
  #} end def

  def ClearGeneSet(self, gene_set_id): #{
    if (None == getattr(self, gene_set_id)): #{
      if ("breakpoint" in gene_set_id): #{
        depth = 4
        list_end = True
      elif ("nearby" in gene_set_id):
        depth = 3
        list_end = False
      else:
        depth = 2
        list_end = True
      #} end if
      setattr(self, gene_set_id, MultiDictCls("none", depth,
        list_end=list_end))
    # for backwards compatibility
    elif (isinstance(getattr(self, gene_set_id), list)):
      setattr(self, gene_set_id, list())
    elif (isinstance(getattr(self, gene_set_id), MultiDictCls)): #{
      getattr(self, gene_set_id).Initialize()
    else:
      raise GroupedCandidateError("Unrecognized gene set class type: %s" %
        type(getattr(self, gene_set_id)))
    #} end if
  #} end def

  def AddGenes(self, region_id, new_genes): #{
    # for backwards compatibility
    if (isinstance(new_genes, set)): #{
      old_genes = self.GetGeneList(region_id)
      if (None != old_genes): #{
        new_genes.update(old_genes)
      #} end if
      self.SetGeneList(region_id, sorted(new_genes))
    else:
      if ("breakpoint" in region_id): #{
        self.AddBreakPointGenes(region_id, new_genes)
      else:
        raise GroupedCandidateError("add non-list main genes "
          "not implemented!")
      #} end if
    #} end if
  #} end def

  def AddBreakPointGenes(self, region_id, new_genes): #{
    #ErrMsg("In AddBreakPointGenes")
    main_region = region_id.replace("breakpoint_","")
    main_genes = self.GetGeneList(main_region)
    #ErrMsg("Main: %s" % main_genes.OutputString(sort=True))
    for id_dict in new_genes: #{
      #ErrMsg("New gene: %s (%s), gene in main: %s, feature in main: %s" %
      #  (id_dict['gene'], id_dict['feature'], id_dict['gene'] in main_genes,
      #  id_dict['feature'] in main_genes[id_dict['gene']]))
      add_transcript = False
      if (id_dict['gene'] in main_genes and
          id_dict['feature'] in main_genes[id_dict['gene']]): #{
        add_transcript = True
      elif (id_dict['transcript'].lower().endswith("non_coding")):
        gene_id = "%s.non_coding" % id_dict['gene']
        #ErrMsg("  Checking non-coding: %s" % gene_id)
        if (gene_id in main_genes and
            id_dict['feature'] in main_genes[gene_id]): #{
          add_transcript = True
        #} end if
      #} end if
      if (add_transcript): #{
        #ErrMsg("  Adding %s.%s feature to candidate!" %
        #  (id_dict['feature'], id_dict['num']))
        keys = [id_dict['gene'], id_dict['transcript'], id_dict['feature']]
        old_gene_list = self.GetGeneList(region_id)
        old_gene_list.Add(keys, id_dict['num'])
      #} end if
    #} end for
  #} end if

  def GetGeneList(self, region_id): #{
    if ("A" == region_id): #{
      return self.genes_A
    elif("nearby_A" == region_id): #{
      return self.nearby_A
    elif ("breakpoint_A" == region_id): #{
      return self.breakpoint_genes_A
    elif ("B" == region_id): #{
      return self.genes_B
    elif("nearby_B" == region_id): #{
      return self.nearby_B
    elif ("breakpoint_B" == region_id): #{
      return self.breakpoint_genes_B
    else:
      raise GroupedCandidateError("invalid region id: %s" % region_id)
    #} end if
  #} end def

  def SetGeneList(self, region_id, new_gene_list): #{
    if ("A" == region_id): #{
      self.genes_A = new_gene_list
    elif ("breakpoint_A" == region_id): #{
      self.breakpoint_genes_A = new_gene_list
    elif ("B" == region_id): #{
      self.genes_B = new_gene_list
    elif ("breakpoint_B" == region_id): #{
      self.breakpoint_genes_B = new_gene_list
    else:
      raise GroupedCandidateError("invalid region id: %s" % region_id)
    #} end if
    if ("breakpoint" in region_id): #{
      self.break_genes = (self.breakpoint_genes_A, self.breakpoint_genes_B)
    elif ("nearby" in region_id):
      self.nearby_genes = (self.nearby_A, self.nearby_B)
    else:
      self.main_genes = (self.genes_A, self.genes_B)
    #} end if
  #} end def

  def OverlapsGene(self, region_id, include_introns=False,
      exclude_UTRs=False, exclude_non_coding=False):
    gene_list = self.GetGeneList(region_id)
    if (isinstance(gene_list, MultiDictCls)): #{
      if ("nearby" in region_id): #{
        include_introns = True
        exclude_UTRs = False
      #} end if
      if (0 == len(gene_list)): #{
        return False
      #} end if
      if (include_introns and not exclude_non_coding and not exclude_UTRs): #{
        return True
      #} end if
      # ensure that the region overlaps a coding and/or non-UTR
      #  gene region
      for gene_name in gene_list: #{
        parts = set()
        if (exclude_non_coding and IsNonCoding(gene_name)): #{
          continue
        elif (include_introns and not exclude_UTRs):
          return True
        elif ("breakpoint" in region_id):
          for transcript in gene_list[gene_name]: #{
            if (exclude_non_coding and IsNonCoding(transcript)): #{
              continue
            #} end if
            parts.update(gene_list[gene_name][transcript])
          #} end for
        else:
          parts.update(gene_list[gene_name])
        #} end if
        if (exclude_UTRs): #{
          parts = filter(IsNotUTR, parts)
        #} end if
        if (not include_introns and INTRON in parts): #{
          parts.remove(INTRON)
        #} end if
        if (0 < len(parts)): #{
          return True
        #} end if
      #} end for
      return False
      # main genes have gene_name(parts_list)
      # nearby genes have gene_name(transcript(distance))
      # breakpoint genes have gene_name(transcript(part(num_list)))
    else:
      if (exclude_UTRs): #{
        gene_list = filter(IsNotUTR, gene_list)
      #} end if
      if (exclude_non_coding): #{
        gene_list = filter(IsNotNonCoding, gene_list)
      #} end if
      if (IsEmpty(gene_list)): #{
        return False
      #} end if
      return (gene_list[0].lower() not in ["none", "no_exons"])
    #} end if
  #} end def

  def GeneSetsOverlap(self): #{
    gene_set_A = set(map(RemoveExonTypeLabel, self.genes_A))
    #print "Before nearby A: %s" % ",".join(gene_set_A)
    gene_set_A.update(map(RemoveExonTypeLabel, self.nearby_A))
    #print "After nearby A: %s" % ",".join(gene_set_A)
    gene_set_B = set(map(RemoveExonTypeLabel, self.genes_B))
    #print "Before nearby B: %s" % ",".join(gene_set_B)
    gene_set_B.update(map(RemoveExonTypeLabel, self.nearby_B))
    #print "After nearby B: %s" % ",".join(gene_set_B)
    in_both = gene_set_A.intersection(gene_set_B)
    #print "Before remove none: %s" % ",".join(in_both)
    in_both.discard("none")
    #print "After remove none: %s" % ",".join(in_both)
    return (0 < len(in_both))
  #} end def

  def BioTypeEventCoordsStr(self, lib, bio_type): #{
    if ("fusion" == bio_type.lower()): #{
      data_A = list([
        self.align_info_A.chrom,
        "%i" % min(self.align_info_A.genome_start,
                   self.align_info_A.genome_end),
        "%i" % max(self.align_info_A.genome_start,
                   self.align_info_A.genome_end),
        lib,
        "%i%s(A)" % (self.group_id, self.candidate_id),
        ",".join(self.genes_A)
      ])
      str_A = "\t".join(data_A)
      data_B = list([
        self.align_info_B.chrom,
        "%i" % min(self.align_info_B.genome_start,
                   self.align_info_B.genome_end),
        "%i" % max(self.align_info_B.genome_start,
                   self.align_info_B.genome_end),
        lib,
        "%i%s(B)" % (self.group_id, self.candidate_id),
        ",".join(self.genes_B)
      ])
      str_B = "\t".join(data_B)
      return "\n".join([str_A, str_B])
    else:
      if ("itd" == bio_type.lower() and self.gap): #{
        data = list([
          self.align_info_B.chrom,
          "%i" % self.align_info_B.genome_start,
          "%i" % self.align_info_B.genome_end,
          lib,
          "%i%s" % (self.group_id, self.candidate_id),
          ",".join(self.genes_A)
        ])
      elif ("itd" == bio_type.lower() or "ptd" == bio_type.lower()): #{
        data = list([
          self.align_info_A.chrom,
          "%i" % min(self.align_info_A.genome_end,
                     self.align_info_B.genome_start),
          "%i" % max(self.align_info_A.genome_end,
                     self.align_info_B.genome_start),
          lib,
          "%i%s" % (self.group_id, self.candidate_id),
          ",".join(self.genes_A)
        ])
      else:
        raise GroupedCandidateError("unrecognized biological type: "
          "\"%s\"" % bio_type)
      #} end if
      return "\t".join(data)
    #} end if
  #} end def

  def ClearRepeats(self): #{
    self.repeats_A = list()
    self.repeats_B = list()
  #} end def

  def AddRepeats(self, repeats_dict, keep_current=False): #{
    if (0 == len(repeats_dict.keys())): #{
      return
    #} end if
    if ("A" in repeats_dict): #{
      new_repeats = set(repeats_dict["A"].values())
      if (keep_current and not IsEmpty(self.repeats_A)): #{
        new_repeats.update(self.repeats_A)
      #} end if
      self.repeats_A = sorted(new_repeats)
    #} end if
    if ("B" in repeats_dict): #{
      new_repeats = set(repeats_dict["B"].values())
      if (keep_current and not IsEmpty(self.repeats_B)): #{
        new_repeats.update(self.repeats_B)
      #} end if
      self.repeats_B = sorted(new_repeats)
    #} end if
  #} end def

  def OverlapsRepeat(self): #{
    #if (1 == len(self.repeats_A)): #{
    #  if ("N/A"  != self.repeats_A[0] and
    #      "none" != self.repeats_A[0].lower()):
    #    return True
    #  #} end if
    #elif (1 < len(self.repeats_A)): #{
    #  return True
    ##} end if
    #if (1 == len(self.repeats_B)): #{
    #  if ("N/A"  != self.repeats_B[0] and
    #      "none" != self.repeats_B[0].lower()):
    #    return True
    #  #} end if
    #elif (1 < len(self.repeats_B)): #{
    #  return True
    #} end if
    if (IsEmpty(self.repeats_A) and IsEmpty(self.repeats_B)): #{
      return False
    #} end if
    return True
  #} end def
#} end class

def ParseFullEventID(object, full_id): #{
  id_pattern = (r"^(?P<group_id>[0-9]+)(?P<candidate_id>[a-z]+)\)?"
    r"(?P<region>[AB])?$")
  id_match = re.search(id_pattern, full_id)
  if (None == id_match): #{
    raise GroupedCandidateError("cannot get member ID from %s" % full_id)
  #} end if
  #object.group_id = int(id_match.group('group_id'))
  #object.candidate_id = id_match.group('candidate_id')
  #if (None != id_match.group("region")): #{
  #  object.region_id = id_match.group("region")
  #} end if
  id_dict = dict()
  id_dict['group_id'] = int(id_match.group('group_id'))
  id_dict['candidate_id'] = id_match.group('candidate_id')
  if (None != id_match.group("region")): #{
    id_dict['region_id'] = id_match.group("region")
  #} end if
  if (None != object): #{
    for id_type in id_dict: #{
      setattr(object, id_type, id_dict[id_type])
    #} end for
  #} end if
  return id_dict
#} end def

def ParseRegionCoords(coords_str, coords_type, warnings, ordered=True,
    allow_na=False, allow_single=False): #{
  rule = "must be positive integer"
  if (allow_na): #{
    rule += " or \"N/A\""
  #} end if
  (start, end) = coords_str.split("-")
  if (allow_na and "N/A-N/A" == coords_str): #{
    return (start, end)
  #} end if
  if (not start.isdigit()): #{
    if (not allow_na or "N/A" != start): #{
      warnings.append(ParseWarningCls("invalid %s start coordinate: " %
        coords_type, start, rule))
    #} end if
  else:
    start = int(start)
  #} end if
  if (not end.isdigit()): #{
    if (not allow_na or "N/A" != start): #{
      warnings.append(ParseWarningCls("invalid %s end coordinate: " %
        coords_type, end, rule))
    #} end if
  else:
    end = int(end)
  #} end if
  if (allow_single): #{
    order_test = (start > end)
  else:
    order_test = (start >= end)
  #} end if
  if (ordered and order_test): #{
    warnings.append(ParseWarningCls("incorrectly ordered %s coordinates" %
      coords_type, coords_str, "start should be less than end"))
  #} end if
  return (start, end)
#} end def

class ContigInfoCls: #{
  def __init__(self, ctg_info_str): #{
    self.ParseCtgInfo(ctg_info_str)
  #} end def

  def ParseCtgInfo(self, ctg_info_str): #{
    ctg_info_pattern = r"(?P<id>.+)\((?P<length>\d+)bp:(?P<k_form>2k[-+]\d+)\)"
    ctg_info_match = re.search(ctg_info_pattern, ctg_info_str)
    self.id     = ctg_info_match.group('id')
    self.length = int(ctg_info_match.group('length'))
    self.k_form = ctg_info_match.group('k_form')
  #} end def

  def ToString(self, with_kform=True): #{
    data_str = "%s(%ibp" % (self.id, self.length)
    if (with_kform): #{
      data_str += ":%s" % self.k_form
    #} end if
    data_str += ")"
    return data_str
  #} end def
#} end class

class AlignmentInfoCls: #{
  def __init__(self, align_info_str): #{
    self.ParseAlignmentInfo(align_info_str)
  #} end def

  def ParseAlignmentInfo(self, align_info_str): #{
    (coords, metrics) = align_info_str.split(";")
    (ctg_coords, genome_coords) = coords.split("=")
    coord_pattern = r"(?P<id>.+):(?P<start>\d+)-(?P<end>\d+)\("
    coord_match = re.search(coord_pattern, ctg_coords)
    self.ctg_start = int(coord_match.group('start'))
    self.ctg_end   = int(coord_match.group('end'))
    self.ctg_coords = CoordPairCls(self.ctg_start, self.ctg_end)
    coord_match = re.search(coord_pattern, genome_coords)
    self.chrom = coord_match.group('id').replace('chr', '').upper()
    self.genome_start = int(coord_match.group('start'))
    self.genome_end   = int(coord_match.group('end'))
    self.gen_coords = CoordPairCls(self.genome_start, self.genome_end)
    (align_fraction, identity) = metrics.split(",")
    self.align_fraction = float(GetFieldValue(align_fraction))
    self.identity       = float(GetFieldValue(identity))
  #} end def

  def ToString(self): #{
    return ("ctg:%s(%ibp)=chr%s:%i-%i(%ibp,%s);AF:%.2f,PID:%.2f" %
        (self.ContigCoordsString(), self.ContigSpan(), self.chrom,
         self.genome_start, self.genome_end, self.GenomeSpan(),
         self.Strand(), self.align_fraction, self.identity))
  #} end def

  def ContigCoordsString(self): #{
    return "%i-%i" % (self.ctg_start, self.ctg_end)
  #} end def

  def ContigSpan(self): #{
    return abs(self.ctg_start - self.ctg_end) + 1
  #} end def

  def GenomeSpan(self): #{
    return abs(self.genome_start - self.genome_end) + 1
  #} end def

  def Strand(self): #{
    if (self.genome_start < self.genome_end): #{
      return "+"
    else:
      return "-"
    #} end if
  #} end def

  def IsPosStrand(self): #{
    return ("+" == self.Strand())
  #} end def

  def IsNegStrand(self): #{
    return not self.IsPosStrand()
  #} end def
#} end class

def IsLeftBreakpoint(is_ctg_start, is_pos_strand): #{
  if (is_ctg_start): #{
    if (is_pos_strand): #{
      return False
    else:
      return True
    #} end if
  else:
    if (is_pos_strand): #{
      return True
    else:
      return False
    #} end if
  #} end if
#} end def

def DetermineBreakpointSide(ctg_part, strand): #{
  if (START == ctg_part): #{
    if ("+" == strand): #{
      return RIGHT
    elif ("-" == strand):
      return LEFT
    else:
      raise GroupedCandidateError("unrecognized alignment strand: \"%s\"" %
        strand)
    #} end if
  elif (END == ctg_part):
    if ("+" == strand): #{
      return LEFT
    elif ("-" == strand):
      return RIGHT
    else:
      raise GroupedCandidateError("unrecognized alignment strand: \"%s\"" %
        strand)
    #} end if
  else:
    raise GroupedCandidateError("unrecognized contig part: \"%s\"" %
      ctg_part)
  #} end if
#} end def

#### EXCEPTION CLASSES ####
class GroupedCandidateError(MyError): #{
  pass
#} end class
