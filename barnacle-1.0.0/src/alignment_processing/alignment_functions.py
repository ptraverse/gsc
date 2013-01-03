#! /usr/bin/env python
"""
alignment_functions.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import time

# import custom modules
from utils.error import MyError
from utils.general import TimeSpent, AddChr, NormalizeChrID
from utils.messages import ErrMsg, LogMsg, DebugMsg#, ExtremeDebugMsg
from parsers.psl_parser import parse as ParsePSL

def ParseAlignmentFile(aligns_path, filters, log_info,
    ParseFunction=ParsePSL): #{
  start = time.time()
  LogMsg(log_info, "Parsing alignment file: %s" % aligns_path)
  aligns = None
  # parse the alignment file
  aligns = ParseFunction(aligns_path, filters, log_info=log_info)
  DebugMsg(log_info, "Found %i alignments" % len(aligns))
  LogMsg(log_info, "Time spent parsing alignment file: %s" % TimeSpent(start))
  if (None == aligns or 0 == len(aligns)): #{
    raise NoAlignmentsError("No alignments were found.")
  #} end if
  return aligns
#} end def

def CalcAlignOverlap(align1, align2): #{
  return CalcOverlap(align1.qstart, align1.qend, align2.qstart, align2.qend)
#} end def

def CalcOverlap(startA, endA, startB, endB): #{
  # order the coordinate pairs so that the smaller coordinate is first
  leftA  = min(startA, endA)
  rightA = max(startA, endA)
  leftB  = min(startB, endB)
  rightB = max(startB, endB)
  # determine the overlap (will be negative if there is a space)
  overlap_left  = max(leftA,  leftB)
  overlap_right = min(rightA, rightB)
  overlap = (overlap_right - overlap_left) + 1
  #ErrMsg("lA: %i, rA: %i, lB: %i, rB: %i, lO: %i, rO: %i, O: %i" %
  #  (leftA, rightA, leftB, rightB, overlap_left, overlap_right, overlap))
  overlap_fraction = 0.0
  # if the coordinates overlapped, calculate the fraction
  if (0 < overlap): #{
    spanA = (rightA - leftA) + 1
    spanB = (rightB - leftB) + 1
    min_span = min(spanA, spanB)
    #ErrMsg("spanA: %i, spanB: %i, min: %i" % (spanA, spanB, min_span))
    overlap_fraction = float(overlap) / float(min_span)
  #} end if
  #ErrMsg("Overlap fraction: %f" % overlap_fraction)
  return (overlap, overlap_fraction)
#} end def

def FixAlign(align, chromosomal=True): #{
  if (hasattr(align, "fixed") and align.fixed): #{
    return align
  #} end if
  if (not hasattr(align, "multi_mapped")): #{
    # assume the alignment was not multi_mapped
    align.multi_mapped = False
  #} end if
  # convert coordinates to integers
  align.qstart   = int(align.qstart)
  align.qend     = int(align.qend)
  align.tstart   = int(align.tstart)
  align.tend     = int(align.tend)
  align.score    = int(align.score)
  align.identity = float(align.identity)
  # assume that the alignment does not contain a query gap
  align.qgap = False
  # assume that the alignment is not spliced
  align.spliced = False
  # ensure that the chromosome name does not contain chr
  if (chromosomal): #{
    align.target = NormalizeChrID(align.target)
  #} end if
  #if (align.target.startswith("chr")): #{
  #  align.target = align.target[3:]
  #} end if
  # ensure that chromosome name contains "chr"
  #if ("chr" != align.target[0:3]): #{
  #  align.target = "chr%s" % align.target
  #} end if
  # convert chrMT to chrM
  #if ("chrMT" == align.target): #{
  #  align.target = "chrM"
  #} end if
  #if ("MT" == align.target): #{
  #  align.target = "M"
  #} end if
  # fix blat alignment
  if ("blat" == align.method): #{
    # convert numeric values to integers
    #align.match      = int(align.match)
    #align.mismatch   = int(align.mismatch)
    #align.query_len  = int(align.query_len)
    #align.target_len = int(align.target_len)
    #align.qnuminsert = int(align.qnuminsert)
    #align.tnuminsert = int(align.tnuminsert)
    numeric_fields = ["match", "mismatch", "query_len", "target_len",
      "qnuminsert", "qbaseinsert", "tnuminsert", "tbaseinsert"]
    for field in numeric_fields: #{
      try:
        int_value = int(getattr(align, field))
        setattr(align, field, int_value)
      except ValueError, e:
        raise AlignmentError("%s must be a numeric value, not \"%s\"" %
          (field, getattr(align, field)))
      #} end try
    #} end for
    # reorder target co-ordinates of alignment to the "-" strand
    if ("-" == align.query_strand): #{
      temp_tstart  = align.tstart
      align.tstart = align.tend
      align.tend   = temp_tstart
    #} end if
    # fix the psl string
    if (None != align.psl_str): #{
      PSL_TNAME_COL = 13
      psl_info = align.psl_str.split("\t")
      target_name = psl_info[PSL_TNAME_COL]
      # convert chrMT to chrM
      # ensure that the psl string has "chr" in the chromosome name
      #if ("chrMT" == target_name): #{
      #  target_name = "chrM"
      #} end if
      #if ("chr" != target_name[0:3]): #{
      #  target_name = "chr%s" % target_name
      #} end if
      if (chromosomal): #{
        target_name = NormalizeChrID(target_name, use_chr=True)
      #} end if
      psl_info[PSL_TNAME_COL] = target_name
      align.psl_str = "\t".join(psl_info)
    #} end if
    # check whether the alignment contains a query gap
    if (0 < align.qnuminsert): #{
      #print "QGap: %i\n" % int(align.qnuminsert)
      #print "PSL: %s\n" % align.psl_str
      align.qgap = True
    #} end if
    # use target gap value to determine whether the alignment is spliced
    if (0 < align.tnuminsert): #{
      align.spliced = True
    #} end if
  # fix exonerate alignments
  elif ("exonerate" == align.method or "exon" == align.method): #{
    # set the query strand
    if (align.target_strand): #{
      align.query_strand = "+"
    else:
      align.query_strand = "-"
      align.blocks.reverse()
    #} end if
    # use the number of blocks to determine whether the alignment is spliced
    if (1 < len(align.blocks)): #{
      align.spliced = True
    #} end if
  #} end if
  # alias for target blocks
  align.target_blocks = align.blocks
  align.fixed = True
  return align
#} end def

def FixNonChromAlign(align): #{
  return FixAlign(align, chromosomal=False)
#} end def

def IntifyBlock(block): #{
  return map(int, block)
#} end def

def WriteBlockCoords(align, id, file, use_chr=False): #{
  target = NormalizeChrID(align.target, use_chr=use_chr)
  #if (use_chr): #{
    #if (align.target.startswith("chr")): #{
    #  target = align.target
    #else:
    #  target = "chr%s" % align.target
    #} end if
  #  target = AddChr(target)
  #else:
    #target = align.target.replace("chr", "", 1)
  #} end if
  for block in map(IntifyBlock, align.target_blocks): #{
    (block_start, block_end) = block
    block_left  = min(block_start, block_end)
    block_right = max(block_start, block_end)
    coord_str = "%s %i %i %i\n" % (target, block_left, block_right, id)
    file.Write(coord_str)
  #} end for
#} end def

#def PrintAligns(self, aligns): #{
#  LogMsg(self, "---")
#  for i in range(len(aligns)): #{
#    LogMsg(self, "%i) " % i, newline=False)
#    PrintAlign(self, aligns[i])
#  #} end for
#  LogMsg(self, "---")
#} end def

def AlignListString(align_list): #{
  align_strs = list()
  align_strs.append("---")
  for i in range(len(align_list)): #{
    align_strs.append("%i) %s" % (i,AlignString(align_list[i])))
  #} end for
  align_strs.append("---")
  return "\n".join(align_strs)
#} end def

#def PrintAlign(self, align): #{
#  ExtremeDebugMsg(self,
#    "S:%i ID:%.3f ML:%i QL:%i Q:%s QS:%i QE:%i T:%s TS:%i TE:%i" %
#    (int(align.score), float(align.identity),
#     int(align.match_len), int(align.query_len),
#     align.query,  align.qstart, align.qend,
#     AddChr(align.target), align.tstart, align.tend))
#} end def

def AlignString(align): #{
  return ("S:%i ID:%.3f ML:%i QL:%i Q:%s QS:%i QE:%i T:%s TS:%i TE:%i" %
    (int(align.score), float(align.identity),
     int(align.match_len), int(align.query_len),
     align.query,  align.qstart, align.qend,
     AddChr(align.target), align.tstart, align.tend))
#} end def

#def PrintShortAlign(self, align): #{
#  ExtremeDebugMsg(self, "S:%i ID:%.1f Q:%s QS:%i QE:%i QL:%i T:%s" %
#    (int(align.score), float(align.identity), align.query, int(align.qstart),
#     int(align.qend), int(align.query_len), AddChr(align.target)))
#} end def

def ShortAlignString(align): #{
  return ("S:%i ID:%.1f Q:%s QS:%i QE:%i QL:%i T:%s" %
    (int(align.score), float(align.identity), align.query, int(align.qstart),
     int(align.qend), int(align.query_len), AddChr(align.target)))
#} end def

class AlignBlockCls: #{
  def __init__(self, block_coords_tuple): #{
    if (block_coords_tuple[0] <= block_coords_tuple[1]): #{
      self.strand = "+"
    else:
      self.strand = "-"
    #} end if
    (self.start, self.end) = sorted(map(int, block_coords_tuple))
    self.span = (self.end-self.start)+1
  #} end def

  def __str__(self): #{
    return "%i-%i" % (self.start, self.end)
  #} end def

  def Gap(self, other): #{
    return min((other.start-self.end),(self.start-other.end))-1
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class AlignmentError(MyError): #{
  pass
#} end class

class NoAlignmentsError(MyError): #{
  pass
#} end class
