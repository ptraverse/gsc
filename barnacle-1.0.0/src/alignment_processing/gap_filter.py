#! /usr/bin/env python
"""
gap_filter.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import re

# import custom modules
from utils.error import MyError
from utils.general import SetupMainClass, NormalizeChrID, NonStandardChr
from utils.messages import LogMsg, DebugMsg, ExtremeDebugMsg
from utils.subprocesses import RunCommandFromString
from gene_overlap import (ChooseBestTranscripts, GenesOverlapped, NearbyGenes)
from alignment_functions import (IntifyBlock, ShortAlignString)

# CONSTANTS
GAP_REALIGNER_FAIL = 0
GAP_REALIGNER_MIN_SUCCESS = 1
GAP_REALIGNER_MAX_SUCCESS = 100
GAP_REALIGNER_ONLY_NON_BASIC = 110

class GapFilterCls: #{
  def __init__(self, ctg_seq_file, gap_out_path, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    # assume no gapped alignments will be found for the current contig
    self.ctg_seq_file = ctg_seq_file
    self.gap_out_path = gap_out_path
    self.ResetFilter()
  #} end def

  def ResetFilter(self): #{
    self.gapped_event_found = "N"
    self.ctg_seq_id = ""
    self.ctg_seq = ""
    self.seq_too_long = False
    self.num_gapped_aligns = 0
    self.max_score = 0
    self.multi_mapped = False
  #} end def

  def FindGappedAlignment(self, curr_align, get_next_seq): #{
    if (self.seq_too_long): #{
      return False
    #} end if
    #LogMsg(self, "Running gap filter") # DEBUG
    if (None == curr_align.psl_str): #{
      msg = "Cannot use gap_realigner on alignment without psl string"
      raise GapFilterError(msg)
    #} end if
    psl_str = curr_align.psl_str
    #DebugMsg(self, psl_str)
    # get the contig info from the contig sequence file
    if (get_next_seq): #{
      self.GetSequenceForContig(curr_align.query)
    #} end if
    #LogMsg(self, "Align: %s" % curr_align.query)
    #LogMsg(self, " Ctgs: %s" % self.ctg_seq_id)
    if (self.ctg_seq_id != curr_align.query): #{
      msg = "Error: contig sequence not in same order as alignments."
      msg += " From sequence file: \"%s\"." % self.ctg_seq_id
      msg += " From alignment: \"%s\"." % curr_align.query
      raise GapFilterError(msg)
    #} end if
    # if the alignment should be passed to the gap realigner
    if (self.FilterCheck(curr_align)): #{
      #DebugMsg(self, "--alignment passed gap filter") # DEBUG
      # run next gap-checking step
      ChooseBestTranscripts(curr_align,
        self.options.transcript_selection_buffer)
      olap_genes = GenesOverlapped(curr_align,
        self.options.add_gene_annotation)
      nearby_genes = NearbyGenes(curr_align, self.options.add_gene_annotation)
      gap_realign_args = [
        self.options.gap_realigner,
        "\"%s\"" % psl_str,
        self.ctg_seq,
        "\"%s NEARBY:%s\"" % (olap_genes, nearby_genes),
        self.gap_out_path,
        self.options.gap_config,
        "%i"   % self.options.min_gap_size,
        "%.2f" % self.options.min_gap_pid,
        "%.2f" % self.options.min_gap_fract,
        "MM:%i" % self.multi_mapped,
      ]
      if (self.options.gap_debug): #{
        gap_realign_args.append("-d")
      #} end if
      gap_realign_cmd = " ".join(gap_realign_args)
      if (self.log_info['debug']): #{
        #LogMsg(self, "   %s" % psl_str) # DEBUG
        ExtremeDebugMsg(self, "Gap realigner Command:\n    %s" %
          gap_realign_cmd)
      #} end if
      if (self.options.arg_max < len(gap_realign_cmd)): #{
        ExtremeDebugMsg(self, "WARNING: Cannot use gap realigner on contig "
          "%s: sequence is too long" % self.ctg_seq_id)
        self.seq_too_long = True
        return False
      #} end if
      #LogMsg(self, "Checking gapped alignment...") # DEBUG
      realigner_return = RunCommandFromString(gap_realign_cmd,
        dpt=self.options.dpt)
      #DebugMsg(self, "Gap realigner returned: %i" % realigner_return) # DEBUG
      if (GAP_REALIGNER_MIN_SUCCESS <= realigner_return and
          GAP_REALIGNER_MAX_SUCCESS >= realigner_return): #{
        if (GAP_REALIGNER_MAX_SUCCESS == realigner_return): #{
          self.more_than_99 = True
        #} end if
        self.num_gapped_aligns += realigner_return
        DebugMsg(self, "%i gapped alignment(s) found for %s" %
          (realigner_return, self.ctg_seq_id))
        self.gapped_event_found = "Y"
        return True
      elif (GAP_REALIGNER_ONLY_NON_BASIC == realigner_return):
        ExtremeDebugMsg(self, "Warning: gap realigner found only "
          "non-basic events.")
      elif (0 > realigner_return):
        raise GapFilterError("Gap realigner command terminated by "
          "signal %i:\n%s" % (-realigner_return, gap_realign_cmd))
      elif (GAP_REALIGNER_FAIL != realigner_return):
        raise GapFilterError("Error running gap realigner: %i:\n%s" %
          (realigner_return, gap_realign_cmd))
      #} end if
    #elif (GAP_FILTER_FAIL != filter_return): #{
    #  msg = "Error running gap filter: %i" % filter_return
    #  raise GapFilterError(msg)
    #else: # DEBUG
    #  LogMsg(self, "--alignment failed gap filter") # DEBUG
    #} end if
    return False
  #} end def

  def GetSequenceForContig(self, ctg_align_id): #{
    #DebugMsg(self, "Getting next conting sequence...")
    # find the contig info with the correct contig id
    ctg_info = self.ctg_seq_file.next()
    if (">" != ctg_info[0]): #{
      msg = "Error reading contig sequence file: %s" % ctg_info
      raise GapFilterError(msg)
    #} end if
    self.ctg_seq_id = ctg_info[1:].split(None,1)[0]
    #DebugMsg(self, "Ctg id from seq file:   %s" % self.ctg_seq_id)
    while (self.ctg_seq_id != ctg_align_id): #{
      LogMsg(self, "WARNING: contig %s not present in alignment file" %
        self.ctg_seq_id)
      # skip the sequence
      self.ctg_seq_file.next()
      # get the next contig id
      try:
        ctg_info = self.ctg_seq_file.next()
      except StopIteration, e:
        msg  = "Error: end of contig sequence file reached before "
        msg += "end of alignment file."
        raise GapFilterError(msg)
      # end try
      if (">" != ctg_info[0]): #{
        msg = "Error reading contig sequence file: %s" % ctg_info
        raise GapFilterError(msg)
      #} end if
      self.ctg_seq_id = ctg_info[1:].split(None,1)[0]
    #} end while
    # get the contig sequence
    self.ctg_seq = self.ctg_seq_file.next()
    self.ctg_seq = self.ctg_seq.rstrip()
  #} end def

  def FilterCheck(self, align): #{
    ExtremeDebugMsg(self, ShortAlignString(align))
    # check that the contig span is wide enough
    if (not self.CheckContigSpan(align)): #{
      return False
    # check to see whether there are enough unaligned bases for a
    # sufficiently sized gap to be possible
    elif (not self.CheckUnalignedBases(align)): #{
      return False
    # check against the max score seen for the current contig
    elif (not self.CheckScoreFraction(align)): #{
      return False
    # check to see whether the alignment has a high enough percent identity
    elif (not self.CheckPercentIdentity(align)): #{
      return False
    # do not check alignments to mitochondrial DNA if the discard
    # mitochondrial DNA flag is set
    elif (not self.CheckMitochondrial(align)): #{
      return False
    # check that the alignment target is a standard chromosome
    # (not random, hap, etc.)
    elif (not self.CheckChromosome(align)): #{
      return False
    # filter by blocks
    return self.CheckBlocks(align)
  #} end def

  def CheckContigSpan(self, align): #{
    contig_span = abs(align.qstart - align.qend) + 1
    cspan_fract = float(contig_span) / float(align.query_len)
    if (self.options.min_gap_check_len_fract > cspan_fract): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: too short contig span: "
        "%i (%.2f)" % (contig_span, cspan_fract))
      return False
    #} end if
    return True
  #} end def

  def CheckUnalignedBases(self, align): #{
    num_unaligned = align.query_len - align.match - align.mismatch
    if (self.options.min_gap_size > num_unaligned): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: too few unaligned bases")
      return False
    #} end if
    return True
  #} end def

  def CheckScoreFraction(self, align): #{
    if (self.max_score < align.score): #{
      self.max_score = align.score
    else:
      score_fract = float(align.score) / float(self.max_score)
      if (self.options.min_score_fract > score_fract): #{
        ExtremeDebugMsg(self, "NO GAP CHECK: score less than max")
        return False
      #} end if
    #} end if
    # do not check alignments with not enough of the contig involved
    # in the alignment (based on score:length ratio)
    score_fract = float(align.score) / float(align.query_len)
    if (self.options.min_gap_check_score_fract > score_fract): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: too little contig aligned")
      return False
    #} end if
    return True
  #} end def

  def CheckPercentIdentity(self, align): #{
    if (self.options.min_gap_check_pid >= align.identity): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: pid too low")
      return False
    #} end if
    return True
  #} end def

  def CheckMitochondrial(self, align): #{
    if (self.options.no_mito and
        #("chrM" == align.target or "M" == align.target)): #{}
        ("M" == NormalizeChrID(align.target))): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: mitochondrial")
      return False
    #} end if
    return True
  #} end def

  def CheckChromosome(self, align): #{
    #if (None == re.search(r"^(chr)?(\d+|[XY]|MT?)$", align.target)): #{}
    if (NonStandardChr(align.target)): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: non-standard chromosome (%s)" %
        align.target)
      return False
    #} end if
    return True
  #} end def

  def CheckBlocks(self, align): #{
    # assume there are no gaps worth looking at
    self.good_gap_found = False
    # get the blocks, as integers, in the correct order
    blocks = map(IntifyBlock, align.query_blocks)
    if ("-" == align.query_strand): #{
      blocks.reverse()
    #} end if
    ExtremeDebugMsg(self, "MIN: %i\n" % self.options.min_gap_size +
      "Num Blocks: %i Strand: %s " % (len(blocks), align.query_strand) +
      "Blocks: %s" % ",".join([BlockString(block) for block in blocks]))
    #} end if
    gap_start = 1
    for block in blocks: #{
      if ("-" == align.query_strand): #{
        block.reverse()
      #} end if
      gap_end = block[0] - 1
      ExtremeDebugMsg(self, "Gap: %i-%i Block: %i-%i" %
        (gap_start, gap_end, block[0], block[1]))
      # check the block order! the block should start after the previous block
      # and should start before it ends
      if (block[0] < gap_start or block[1] < block[0]): #{
        ExtremeDebugMsg(self,
          "gap_start: %i, block_start: %i, block_end: %i" %
          (gap_start, block[0], block[1]))
        raise GapFilterError("incorrect block ordering: %s" %
          ",".join([BlockString(block) for block in blocks]))
      #} end if
      if (not self.CheckGap(gap_start, gap_end, align.query_len)): #{
        return False
      #} end if
      gap_start = block[1] + 1
    #} end for
    # look for gap at the very end of the contig
    ExtremeDebugMsg(self, "Last Gap: %i-%i" % (gap_start, align.query_len))
    if (not self.CheckGap(gap_start, align.query_len, align.query_len)): #{
      return False
    #} end if
    if (not self.good_gap_found): #{
      ExtremeDebugMsg(self, "NO GAP CHECK: no good gap found")
    #} end if
    return self.good_gap_found
  #} end def

  def OldCheckBlocks(self, align): #{
    self.good_gap_found = False
    #block_sizes = [0]
    #block_sizes.extend(align.block_sizes.rstrip(",").split(","))
    #block_sizes.append(0)
    #query_starts = [0]
    #query_starts.extend(align.qstarts.rstrip(",").split(","))
    #query_starts.append(align.query_len)
    #(prev_size, prev_start) = (None, None)
    blocks = map(IntifyBlock, align.query_blocks)
    if ("-" == align.query_strand): #{
      blocks.reverse()
    #} end if
    ExtremeDebugMsg(self, "\n".join(["MIN: %i" % self.options.min_gap_size,
      "Num Blocks: %i Strand: %s" % (len(blocks), align.query_strand),
      "Blocks: %s" % ",".join(["%i-%i" %
        tuple(sorted([block[0],block[1]])) for block in blocks])]))
    gap_start = 1
    #for bsize, qstart in zip(block_sizes, query_starts): #{
    for block in blocks: #{
      if ("-" == align.query_strand): #{
        block.reverse()
      #} end if
      gap_end = block[0] - 1
      ExtremeDebugMsg(self, "Gap: %i-%i Block: %i-%i" %
        (gap_start, gap_end, block[0], block[1]))
      # check the block order! the block should start after the previous block
      # and should start before it ends
      if (block[0] < gap_start or block[1] < block[0]): #{
        ExtremeDebugMsg(self,
          "gap_start: %i, block_start: %i, block_end: %i" %
          (gap_start, block[0], block[1]))
        raise GapFilterError("incorrect block ordering: %s" %
          ",".join(["%i-%i" % tuple(sorted([block[0],block[1]])) for
          block in blocks]))
      #} end if
      #if (None != prev_size and None != prev_start): #{
      gapsize = (gap_end - gap_start) + 1
      #DebugMsg(self, "PS: %i PL: %i PE: %i CS: %s\n" %
      #  (prev_start, prev_size, prev_start+prev_size, qstart) +
      #  "GAPSIZE: %i" % gapsize)
      if (gapsize >= self.options.min_gap_size): #{
        gap_seq = self.ctg_seq[gap_start-1:gap_end].lower()
        ExtremeDebugMsg(self, "Gap Start: %s Gap End: %s\nSequence: %s" %
          (gap_start, gap_end, gap_seq))
        # make sure that the gap does not include a single 'n'
        pattern = r"(^n[^n]\|[^n]n[^n]\|[^n]n$)"
        if (None != re.search(pattern, gap_seq)): #{
          return False
        #} end if
        # make sure the gap is not just T's or A's at the beginning
        # or end of the contig
        if ((1 < gap_start and gap_end < align.query_len) or
            (None == re.search(r"^a+$", gap_seq) and
             None == re.search(r"^t+$", gap_seq))):
          self.good_gap_found = True
        else:
          ExtremeDebugMsg(self, "Poly-A tail")
      #} end if
      gap_start = block[1] + 1
    #} end for
    # look for gap at the very end of the contig
    gapsize = (int(align.query_len) - gap_start) + 1
    if (gapsize >= self.options.min_gap_size): #{
      gap_seq = self.ctg_seq[gap_start-1:int(align.query_len)].lower()
      ExtremeDebugMsg(self, "Gap Start: %s Gap End: %s\nSequence: %s" %
        (gap_start, align.query_len, gap_seq))
      # make sure the gap is not just T's or A's at the beginning
      # or end of the contig
      if (None == re.search(r"^a+$", gap_seq) and
          None == re.search(r"^t+$", gap_seq)):
        self.good_gap_found = True
      else:
        ExtremeDebugMsg(self, "Poly-A tail")
    #} end if
    #if (self.log_info['debug']): #{
    #  LogMsg(self, "BLOCK SIZES: %s" % align.block_sizes)
    #  for bsize in block_sizes: #{
    #    LogMsg(self, "  %s" % bsize)
    #  #} end for
    #  LogMsg(self, "QUERY STARTS: %s" % align.qstarts)
    #  for qstart in query_starts: #{
    #    LogMsg(self, "  %s" % qstart)
    #  #} end for
    #} end if
    return self.good_gap_found
  #} end def

  def CheckGap(self, gap_start, gap_end, query_len): #{
    gapsize = (gap_end - gap_start) + 1
    if (gapsize >= self.options.min_gap_size): #{
      gap_seq = self.ctg_seq[gap_start-1:gap_end].lower()
      ExtremeDebugMsg(self, "Gap Start: %i Gap End: %i\nSequence: %s" %
        (gap_start, gap_end, gap_seq))
      # make sure that the gap does not include a single 'n'
      pattern = r"(^n[^n]\|[^n]n[^n]\|[^n]n$)"
      if (None != re.search(pattern, gap_seq.lower())): #{
        ExtremeDebugMsg(self, "NO GAP CHECK: single 'n' in gap sequence")
        return False
      #} end if
      # make sure the gap is not just T's or A's at the beginning
      # or end of the contig
      if ((1 < gap_start and gap_end < query_len) or
          (None == re.search(r"^a+$", gap_seq) and
           None == re.search(r"^t+$", gap_seq))):
        self.good_gap_found = True
      else:
        ExtremeDebugMsg(self, "Poly-A tail")
    else:
      ExtremeDebugMsg(self, "Gap too small: %i, Min: %i" %
        (gapsize, self.options.min_gap_size))
    #} end if
    return True
  #} end def
#} end class

def BlockString(block): #{
  return "%i-%i" % tuple(sorted([block[0],block[1]]))
#} end def

#### EXCEPTION CLASSES ####
class GapFilterError(MyError): #{
  """Exception raised for errors encountered by the GapFilterCls class"""
  pass
#} end class
