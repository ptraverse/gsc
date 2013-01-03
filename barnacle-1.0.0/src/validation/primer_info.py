#! /usr/bin/env python
"""
primer_info.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
from utils.messages import LogMsg, DebugMsg
from alignment_processing.alignment_functions import CalcOverlap

# CONSTANTS

class PrimerInfoCls: #{
  def __init__(self, member, event_type, log_info=None): #{
    self.event_type = event_type
    self.group_id   = member.group_id
    self.ctg_id     = member.contig_info.id
    self.ctg_len    = member.contig_info.length
    self.c2g_aligns = [
      [member.align_info_A.ctg_start, member.align_info_A.ctg_end],
      [member.align_info_B.ctg_start, member.align_info_B.ctg_end]]
    if ("itd" in event_type): #{
      (rqstart, rqend, rtstart, rtend) = member.RealignmentCoords(trim=True)
      self.c2g_aligns[1:] = [[rtstart, rtend], [rqstart,rqend]]
      self.internal_gap = member.GapIsInternal()
    #} end if
    self.c2t_aligns = [None, None]
    self.c2t_aligns_extra = [list(), list()]
    self.log_info = log_info
  #} end def

  def EventID(self): #{
    return "G%i_%s" % (self.group_id, self.event_type)
  #} end def

  def ToString(self): #{
    return "G%i_%s: %s(%ibp) %i-%i; %i-%i" % (self.group_id, self.event_type,
      self.ctg_id, self.ctg_len, self.c2g_aligns[0][0], self.c2g_aligns[0][1],
      self.c2g_aligns[1][0], self.c2g_aligns[1][1])
  #} end def

  def C2GAlignOverlap(self): #{
    return (self.c2g_aligns[0][1] - self.c2g_aligns[1][0]) + 1
  #} end def

  def C2TAlignOverlap(self): #{
    return (self.c2t_aligns[0].ctg_end - self.c2t_aligns[1].ctg_start) + 1
  #} end def

  def UpdateC2TAligns(self, full_align): #{
    new_align = PrimerInfoAlignmentCls(full_align)
    #DebugMsg(self, "ALIGN: %s" % new_align.ToString())
    overlaps = [0,0]
    # determine which contig region this alignment is for
    for region_id in [0,1]: #{
      # CalcOverlap returns (overlap, overlap_fraction)
      overlaps[region_id] = CalcOverlap(self.c2g_aligns[region_id][0],
        self.c2g_aligns[region_id][1], new_align.ctg_start, new_align.ctg_end)
      #DebugMsg(self, "  %i overlap: %i" % (region_id, overlaps[region_id][0]))
    #} end for
    align_region = 0
    # if the region 0 overlap is less than the region 1 overlap
    if (overlaps[0][0] < overlaps[1][0]): #{
      align_region = 1
    #} end if
    if (None == self.c2t_aligns[align_region] or
        self.c2t_aligns[align_region].match < new_align.match): #{
      self.c2t_aligns[align_region] = new_align
      self.c2t_aligns_extra[align_region] = list()
    elif (self.c2t_aligns[align_region].match == new_align.match):
      if (new_align == self.c2t_aligns[align_region]): #{
        DebugMsg(self, "New alignment identical to existing!\n%s\n%s" %
          (self.c2t_aligns[align_region], new_align))
        return
      #} end if
      for extra_align in self.c2t_aligns_extra[align_region]: #{
        if (new_align == extra_align): #{
          DebugMsg(self, "New alignment identical to extra!\n%s\n%s" %
            (extra_align, new_align))
          return
        #} end if
      # end for
      #DebugMsg(self, "CHECK NEW_ALIGN NOT IDENTICAL TO ANY EXISTING!")
      DebugMsg(self, "Adding extra alignment to %s" % new_align.transcript)
      self.c2t_aligns_extra[align_region].append(new_align)
    #} end if
  #} end def

  def ChooseAlignments(self): #{
    if ("fusion" in self.event_type or "itd" in self.event_type): #{
      return
    #} end if
    if ("ptd" in self.event_type): #{
      DebugMsg(self, "Choosing alignments for %s" % self.EventID())
      if (0 == len(self.c2t_aligns_extra[0]) or
          0 == len(self.c2t_aligns_extra[1]) or
          (self.c2t_aligns[0].transcript == self.c2t_aligns[1].transcript and
           self.c2t_aligns[0].strand == self.c2t_aligns[1].strand)): #{
        DebugMsg(self, "Using first choice...")
        return
      #} end if
      for extra_align in self.c2t_aligns_extra[1]: #{
        if (self.c2t_aligns[0].transcript == extra_align.transcript and
            self.c2t_aligns[0].strand == extra_align.strand): #{
          DebugMsg(self, "Swapping alignment 1...")
          self.c2t_aligns[1] = extra_align
          return
        #} end if
      #} end for
      for extra_align in self.c2t_aligns_extra[0]: #{
        if (self.c2t_aligns[1].transcript == extra_align.transcript and
            self.c2t_aligns[1].strand == extra_align.strand): #{
          DebugMsg(self, "Swapping alignment 0...")
          self.c2t_aligns[0] = extra_align
          return
        #} end if
      #} end for
      for extra_align0 in self.c2t_aligns_extra[0]: #{
        for extra_align1 in self.c2t_aligns_extra[1]: #{
          if (extra_align0.transcript == extra_align1.transcript and
              extra_align0.strand == extra_align1.strand): #{
            DebugMsg(self, "Swapping both alignments...")
            self.c2t_aligns[0] = extra_align0
            self.c2t_aligns[1] = extra_align1
            return
          #} end if
        #} end for
      #} end for
    #} end if
    raise MyError("ChooseAlignments not implemented!")
  #} end def
#} end class

class PrimerInfoAlignmentCls: #{
  def __init__(self, align): #{
    self.ctg_id = align.query
    self.transcript = align.target
    self.match = align.match
    self.strand = align.query_strand
    self.ctg_start  = align.qstart
    self.ctg_end    = align.qend
    self.ctg_len    = align.query_len
    self.transcript_start = align.tstart
    self.transcript_end   = align.tend
    self.transcript_len   = align.target_len
    self.num_query_gaps = align.qnuminsert
    # only save the blocks if there are query gaps that might need checking
    self.query_blocks = None
    self.target_blocks = None
    if (0 < self.num_query_gaps): #{
      self.query_blocks = align.query_blocks
      self.target_blocks = align.target_blocks
      # if aligned to the negative strand, need to reverse all the blocks
      if ("-" == align.query_strand): #{
        self.query_blocks.reverse()
        self.target_blocks.reverse()
        #for index in range(len(self.query_blocks)): #{
        #  self.query_blocks[index].reverse()
        #} end for
        map(lambda block: block.reverse(), self.query_blocks)
        map(lambda block: block.reverse(), self.target_blocks)
      #} end if
    #} end if
  #} end def

  def __eq__(self, other): #{
    if (None == other): #{
      return False
    #} end if
    if (self.ctg_id           != other.ctg_id or
        self.transcript       != other.transcript or
        self.match            != other.match or
        self.strand           != other.strand or
        self.ctg_start        != other.ctg_start or
        self.ctg_end          != other.ctg_end or
        self.transcript_start != other.transcript_start or
        self.transcript_end   != other.transcript_end or
        self.num_query_gaps   != other.num_query_gaps): #{
      return False
    #} end if
    if (0 < self.num_query_gaps): #{
      if (len(self.query_blocks)  != len(other.query_blocks) or
          len(self.target_blocks) != len(other.target_blocks)): #{
        return False
      #} end if
      for qindex in range(len(self.query_blocks)): #{
        if (self.query_blocks[qindex][0] != other.query_blocks[qindex][0] or
            self.query_blocks[qindex][1] != other.query_blocks[qindex][1]): #{
          return False
        #} end if
      #} end for
      for tindex in range(len(self.target_blocks)): #{
        if (self.target_blocks[tindex][0] !=
            other.target_blocks[tindex][0] or
            self.target_blocks[tindex][1] !=
            other.target_blocks[tindex][1]): #{
          return False
        #} end if
      #} end for
    #} end if
    return True
  #} end def

  def __ne__(self, other): #{
    return not self.__eq__(other)
  #} end def

  def __str__(self): #{
    return self.ToString()
  #} end def

  def ToString(self): #{
    data_str = "M:%i C:%i-%i T(%s):%i-%i Strand:%s" % (self.match,
      self.ctg_start, self.ctg_end, self.transcript,
      self.transcript_start, self.transcript_end, self.strand)
    if (None != self.query_blocks): #{
      data_str += ": %s" % ",".join(["%i-%i" % (block[0],block[1]) for
        block in self.query_blocks])
    #} end if
    if (None != self.target_blocks): #{
      data_str += "::%s" % ",".join(["%i-%i" % (block[0],block[1]) for
        block in self.target_blocks])
    #} end if
    return data_str
  #} end if
#} end class

class PrimerSeqCls: #{
  def __init__(self, id): #{
    self.id    = id
    self.left  = None
    self.right = None
    self.seq   = None
  #} end def
#} end def
