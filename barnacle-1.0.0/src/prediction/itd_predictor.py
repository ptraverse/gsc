#! /usr/bin/env python
"""
itd_predictor.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# TODO
# deal with always excluding non_coding genes when getting ITDs

# import standard modules

# import custom modules
from utils.error import MyError
from utils.general import (MaxMatchLength, IsHomopolymerSequence,
  ReverseComplement)
from utils.messages import LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import GetFilePath, FileBoxCls
from base_predictor import BasePredictorCls

# constants
MIN_DUP = 4
MIN_FULL_CTG_DUP_FRACTION = 0.25
MAX_LOCAL_INVERSION = 50
BIG_GAP_SPAN = 15

class ITDPredictorCls(BasePredictorCls): #{
  def __init__(self, options, log_info=None): #{
    BasePredictorCls.__init__(self, 'itd', 'itd',
      'internal tandem duplication', options, log_info=log_info)
    self.good_topologies = set([
      'gap-tandem-duplication',
      'gap-nontandem-duplication',
    ])
    if (options.allow_non_gap_itds): #{
      self.good_topologies.update([
        'end-duplication',
        'intrachr-non-colinear',
        'junction-duplication',
      ])
    #} end if
    #self.store_seq = True
    #multi_exon_path = GetFilePath(options.output_dir, options.lib,
    #  "%s.multi_exon" % self.ext)
    #self.multi_exon_file = FileBoxCls(multi_exon_path, "w", "cannot "
    #  "create multiple exon duplication %s output file" % self.description)
    #full_ctg_dup_path = GetFilePath(options.output_dir, options.lib,
    #  "%s.full_ctg_dup" % self.ext)
    #self.full_ctg_dup_file = FileBoxCls(full_ctg_dup_path, "w", "cannot "
    #  "create full contig duplication %s output file" % self.description)
    #edge_gap_path = GetFilePath(options.output_dir, options.lib,
    #  "%s.edge_gap" % self.ext)
    #self.edge_gap_file = FileBoxCls(edge_gap_path, "w", "cannot create "
    #  "edge-gap duplication %s output file" % self.description)
    self.AddOutFile("internal_me", "%s.multi_exon" % self.ext,
      "multiple exon duplication %s" % self.description)
    self.AddOutFile("edge", "%s.edge_gap" % self.ext,
      "edge-gap duplication %s" % self.description)
    self.AddOutFile("edge_me", "%s.edge_gap.multi_exon" % self.ext,
      "edge-gap multiple exon duplication %s" % self.description)
    self.AddOutFile("full", "%s.full_ctg_dup" % self.ext,
      "full contig duplication %s" % self.description)
    self.AddOutFile("full_me", "%s.full_ctg_dup.multi_exon" % self.ext,
      "multiple exon full contig duplication %s" % self.description)
    self.num_over_aligned = 0
  #} end def

  def TestGroup(self, group): #{
    DebugMsg(self, "Group %s. Topologies: %s" %
      (group.id, ",".join(group.topologies)))
    # check alignment topology
    if (0 == len(self.good_topologies.intersection(group.topologies))): #{
      DebugMsg(self, "wrong alignment topology")
      return False
    #} end if
    return True
  #} end def

  def TestMember(self, member): #{
    DebugMsg(self, "  Member %s. Topology: %s" %
      (member.candidate_id, member.topology))
    # ignore members with the wrong topologies
    if (member.topology not in self.good_topologies): #{
      DebugMsg(self, "  %s: Wrong alignment topology: %s" %
        (member.IDString(), member.topology))
      return False
    #} end if
    # if the event matches two exon boundaries, it is not an ITD
    if (2 == member.meta_fields['exon_bounds']): #{
      DebugMsg(self, "  %s: Too many exon boundaries" % member.IDString())
      return False
    #} end if
    DebugMsg(self, "  valid number of exon boundaries")
    # if the member has a gap topology
    if (member.gap): #{
      # if the event coordinates are not internal to the contig
      if (not member.GapIsInternal()): #{
        DebugMsg(self, "  gap is not internal!")
        # then the member might not be an ITD
        if (self.options.require_internal_gaps): #{
          return False
        #} end if
        # check the edge gap fraction
        if (None != self.options.min_edge_gap_fraction): #{
          align_len    = float(member.align_info_B.ContigSpan())
          gap_len      = float(member.GapSpan())
          gap_fraction = align_len / gap_len
          if (1 < gap_fraction): #{
            self.num_over_aligned += 1
            DebugMsg(self, "  WARNING: gapped candidate with aligned length "
               "greater than gap length! (%s)" % member.IDString())
            ExtremeDebugMsg(self, "  %s" % member.DataString())
          #} end if
          # if not enough of the gap is aligned, it is not an ITD
          if (self.options.min_edge_gap_fraction > gap_fraction): #{
            DebugMsg(self, "  not enough gap realigned!")
            return False
          #} end if
          DebugMsg(self, "  sufficient gap realignment")
        #} end if
      #} end if
      DebugMsg(self, "  member with internal gap")
    #} end if
    # if the member overlaps an annotated repeat region,
    # it might not be an ITD
    if (self.options.filter_itd_repeats and (member.OverlapsRepeat())): #{
      DebugMsg(self, "  member overlaps repeat!")
      return False
    #} end if
    # if the gene set is not empty, then the member might be an ITD
    return self.ITDGeneCheck(member)
    #if (member.ITDGeneCheck(exclude_non_coding=True)): #{
    #  self.bio_type_msg += ", %s" % member.bio_type_msg
    #  return True
    #} end if
  #} end def

  def ITDGeneCheck(self, member): #{
    if (member.gap): #{
      if (member.OverlapsGene("breakpoint_A", exclude_UTRs=True,
          exclude_non_coding=self.options.exclude_non_coding)):
        DebugMsg(self, "  overlapping a gene")
        return True
      #} end if
      DebugMsg(self, "  not overlapping any gene")
    else:
      # check that neither gene set is empty
      if (not member.OverlapsGene("breakpoint_A", exclude_UTRs=True,
          exclude_non_coding=self.options.exclude_non_coding) or
          not member.OverlapsGene("breakpoint_B", exclude_UTRs=True,
          exclude_non_coding=self.options.exclude_non_coding)): #{
        DebugMsg(self, "  empty gene set")
        return False
      #} end if
      # check that the gene sets overlap
      if (member.GeneSetsOverlap()): #{
        DebugMsg(self, "  with overlapping gene sets")
        return True
      #} end if
      DebugMsg(self, "  with disjoint gene sets")
    #} end if
    DebugMsg(self, "  failed internal tandem duplication gene check")
    return False
  #} end def

  def OutputEvent(self, event): #{
    multi_exon = False
    # write edge-gap events to separate files
    if (self.write_special): #{
      for member in event.members: #{
        if (1 < len(member.blocks_B)): #{
          DebugMsg(self, "Multiple-exons duplication event member")
          multi_exon = True
        #} end if
        if (member.gap and not member.GapIsInternal()): #{
          if (MIN_FULL_CTG_DUP_FRACTION <=
              min(member.align_info_A.align_fraction,
                  member.align_info_B.align_fraction)): #{
            DebugMsg(self, "Full contig duplication event")
            # write the event to the appropriate full contig duplication
            # output file
            if (multi_exon): #{
              self.outfiles["full_me"].Write(event.FullDataString())
            else:
              self.outfiles["full"].Write(event.FullDataString())
            #} end if
          else:
            DebugMsg(self, "Edge-gap duplication event")
            # write the event to the appropriate edge-gap output file
            if (multi_exon): #{
              self.outfiles["edge_me"].Write(event.FullDataString())
            else:
              self.outfiles["edge"].Write(event.FullDataString())
            #} end if
          #} end if
          return
        #} end if
      #} end for
      if (multi_exon): #{
        DebugMsg(self, "Multiple-exons internal duplication event")
        # write the event to the multiple-exons output file
        self.outfiles["internal_me"].Write(event.FullDataString())
        return
      #} end if
    #} end if
    BasePredictorCls.OutputEvent(self, event)
  #} end def

  #def ReprocessMember(self, member, contig_alignment, contig_seq): #{}
  def ReprocessMember(self, member, realigned_contig): #{
    contig_alignment = realigned_contig.best_align
    ctg_seq = realigned_contig.sequence
    if (None == contig_alignment): #{
      LogMsg(self, "  No alignments for contig %s" % member.contig_info.id)
      return True
    #} end if
    if (member.IDString() in realigned_contig.full_dup_aligns): #{
      full_dup_coords = realigned_contig.full_dup_regions[member.IDString()]
      (full_start,full_end) = (full_dup_coords.min-1, full_dup_coords.max)
      full_dup_seq = ctg_seq[full_start:full_end]
      # check for microrepeat expansion
      dup_coords = member.align_info_B.ctg_coords
      (start,end) = (dup_coords.min-1, dup_coords.max)
      dup_seq = ctg_seq[start:end]
      while ((full_dup_seq + dup_seq) in ctg_seq): #{
        full_dup_seq += dup_seq
      #} end while
      DebugMsg(self, "Full duplicated sequence: %s" % full_dup_seq)
      targets = sorted(realigned_contig.full_dup_aligns[member.IDString()])
      for transcript_id in targets: #{
        t_seq = self.transcript_sequences[transcript_id]
        if (full_dup_seq in t_seq or ReverseComplement(full_dup_seq) in t_seq): #{
          DebugMsg(self, "  False positive: full duplicate alignment (%s)" %
            transcript_id)
          return False
        #} end if
      #} end for
    #} end if
    align_span = contig_alignment.Span()
    DebugMsg(self, "  Realign:%s" % contig_alignment.ToString())
    if (contig_alignment.perfect): #{
      DebugMsg(self, "  False positive: perfect alignment (%i)" %
        contig_alignment.match)
      return False
    #} end if
    if (member.gap and member.GapIsInternal() and
        0 == contig_alignment.num_query_gaps and
        contig_alignment.ctg_start <= member.meta_fields['gap_start'] and
        contig_alignment.ctg_end >= member.meta_fields['gap_end']): #{
      DebugMsg(self, "  False positive: ungapped alignment: %s" %
        contig_alignment.ToString())
      return False
    #} end if
    t_seq = self.transcript_sequences[contig_alignment.target]
    if ("-" == contig_alignment.strand): #{
      t_seq = ReverseComplement(t_seq)
    #} end if
    # do not count "duplications" involving mismatches with the ref
    #if (member.event_seq not in t_seq): #{
    #  LogMsg(self, "Member %s ITD event sequence %s not in transcript %s" %
    #    (member.IDString(), member.event_seq, contig_alignment.target))
    #} end if
    # check any "gap" before the alignment
    DebugMsg(self, "BEFORE alignment:")
    gap_span = (contig_alignment.ctg_start - member.StartOfContig())
    if (MIN_DUP <= gap_span): #{
      if (self.CheckStartAligns(realigned_contig.start_aligns,
          contig_alignment.ctg_start, contig_alignment.target,
          contig_alignment.strand)): #{
        DebugMsg(self, "Local inversion false-positive")
        return False
      #} end if
      gap_start = member.StartOfContig()-1
      gap_end   = contig_alignment.ctg_start-1
      gap_seq   = ctg_seq[gap_start:gap_end]
      buff = (gap_span / BIG_GAP_SPAN)/3
      align_seq = ctg_seq[gap_end:]
      DebugMsg(self, "Contig: %s\n Gap: %s\n Aligned: %s" %
        (ctg_seq, gap_seq, align_seq))
      DebugMsg(self, "Gap_Span:%i, buff:%i" % (gap_span, buff))
      if (BIG_GAP_SPAN < gap_span): #{
        if (gap_seq in align_seq): #{
          DebugMsg(self, "  Good itd prediction: fully duplicated big gap!")
          return True
        #} end if
        if (gap_seq[buff:] in align_seq or gap_seq[:-buff] in align_seq): #{
          DebugMsg(self, "  Good itd prediction: mostly duplicated big gap!\n"
            "  %s (%i)" % (gap_seq, buff))
          return True
        #} end if
      #} end if
      #gap_seq_rc = ReverseComplement(gap_seq)
      #if (gap_seq_rc in ctg_seq or gap_seq_rc in t_seq): #{}
      if (not self.MicroInversionCheck(member, gap_seq, ctg_seq, t_seq) and
          not IsHomopolymerSequence(gap_seq)): #{
        target_start = contig_alignment.ctg_start-1
        target_end   = min(target_start + gap_span,
          realigned_contig.ctg_length)
        target_seq = ctg_seq[target_start:target_end]
        best_match = MaxMatchLength(gap_seq, target_seq)
        min_span = min(align_span, gap_span)
        min_match = max(MIN_DUP,
          int(min_span*self.options.min_edge_gap_fraction))
        DebugMsg(self, "BEFORE\nGap:  %s\nTarg: %s\nMatch: %i Min: %i" %
          (gap_seq, target_seq, best_match, min_match))
        if (min_match <= best_match): #{
          DebugMsg(self, "  Good itd prediction: duplicated gap!")
          return True
        #} end if
      #} end if
    #} end if
    # check any "gap" after the alignment
    gap_span = (member.EndOfContig() - contig_alignment.ctg_end)
    DebugMsg(self, "AFTER alignment:")
    if (MIN_DUP <= gap_span): #{
      if (self.CheckEndAligns(realigned_contig.end_aligns,
          contig_alignment.ctg_end, contig_alignment.target,
          contig_alignment.strand)): #{
        DebugMsg(self, "Local inversion false-positive")
        return False
      #} end if
      gap_start = contig_alignment.ctg_end
      gap_end   = member.EndOfContig()
      gap_seq   = ctg_seq[gap_start:gap_end]
      buff = (gap_span / BIG_GAP_SPAN)/3
      align_seq = ctg_seq[:gap_start]
      DebugMsg(self, "Contig: %s\n Gap: %s\n Aligned: %s" %
        (ctg_seq, gap_seq, align_seq))
      DebugMsg(self, "Gap_Span:%i, buff:%i" % (gap_span, buff))
      if (BIG_GAP_SPAN < gap_span): #{
        if (gap_seq in align_seq): #{
          DebugMsg(self, "  Good itd prediction: fully duplicated big gap!")
          return True
        #} end if
        if (gap_seq[buff:] in align_seq or gap_seq[:-buff] in align_seq): #{
          DebugMsg(self, "  Good itd prediction: mostly duplicated big gap!\n"
            "  %s (%i)" % (gap_seq, buff))
          return True
        #} end if
      #} end if
      #gap_seq_rc = ReverseComplement(gap_seq)
      #if (gap_seq_rc in ctg_seq or gap_seq_rc in t_seq): #{
      if (not self.MicroInversionCheck(member, gap_seq, ctg_seq, t_seq) and
          not IsHomopolymerSequence(gap_seq)): #{
        target_end   = contig_alignment.ctg_end
        target_start = max(0, target_end - gap_span)
        target_seq   = ctg_seq[target_start:target_end]
        best_match = MaxMatchLength(gap_seq, target_seq)
        min_span = min(align_span, gap_span)
        min_match = max(MIN_DUP,
          int(min_span*self.options.min_edge_gap_fraction))
        DebugMsg(self, "AFTER\nGap:  %s\nTarg: %s\nMatch: %i Min: %i" %
          (gap_seq, target_seq, best_match, min_match))
        if (min_match <= best_match): #{
          DebugMsg(self, "  Good itd prediction: duplicated gap!")
          return True
        #} end if
      #} end if
    #} end if
    # check any gaps within the alignment
    if (0 < contig_alignment.num_query_gaps): #{
      prev_right = contig_alignment.query_blocks[0][1]
      for block in contig_alignment.query_blocks[1:]: #{
        DebugMsg(self, "P:%i B:%i-%i" % (prev_right, block[0], block[1]))
        gap_span = (block[0] - prev_right) - 1
        if (MIN_DUP <= gap_span): #{
          DebugMsg(self, "Gap: %i-%i" % (prev_right+1, block[0]-1))
          gap_start = prev_right
          gap_end   = block[0]-1
          if (self.CheckLocalInversion(gap_start, gap_end,
              realigned_contig, contig_alignment.target,
              contig_alignment.strand)): #{
            DebugMsg(self, "Local inversion false-positive")
            return False
          #} end if
          gap_seq   = ctg_seq[gap_start:gap_end]
          #gap_seq_rc = ReverseComplement(gap_seq)
          #if (gap_seq_rc in ctg_seq or
          #    gap_seq_rc in t_seq): #{
          if (self.MicroInversionCheck(member, gap_seq, ctg_seq, t_seq)): #{
            prev_right = block[1]
            continue
          #} end if
          if (IsHomopolymerSequence(gap_seq)): #{
            DebugMsg(self, "  Skipping homopolymer gap: %s" % gap_seq)
            prev_right = block[1]
            continue
          #} end if
          # check before the gap
          target_end   = gap_start
          target_start = max(0, target_end - gap_span)
          target_seq   = ctg_seq[target_start:target_end]
          best_match = MaxMatchLength(gap_seq, target_seq)
          DebugMsg(self, "INTERNAL BEFORE\nGap:  %s\nTarg: %s\nMatch: %i" %
            (gap_seq, target_seq, best_match))
          if (MIN_DUP <= best_match): #{
            DebugMsg(self, "  Good itd prediction: duplicated gap!")
            return True
          #} end if
          # check after the gap
          target_start = block[0]-1
          target_end   = min(target_start + gap_span,
            realigned_contig.ctg_length)
          target_seq   = ctg_seq[target_start:target_end]
          best_match = MaxMatchLength(gap_seq, target_seq)
          DebugMsg(self, "INTERNAL AFTER\nGap:  %s\nTarg: %s\nMatch: %i" %
            (gap_seq, target_seq, best_match))
          if (MIN_DUP <= best_match): #{
            DebugMsg(self, "  Good itd prediction: duplicated gap!")
            return True
          #} end if
        #} end if
        prev_right = block[1]
      #} end for
    #} end if
    # KLUDGE
    if (not member.gap): #{
      LogMsg(self, "  WARNING: automatically passing non-gap ITD prediction: "
        "%s" % realigned_contig.id)
      return True
    #} end if
    DebugMsg(self, "  No gaps have flanking duplications")
    return False
  #} end def

  def MicroInversionCheck(self, member, gap_seq, ctg_seq, t_seq): #{
    dup_span = member.align_info_B.ContigSpan()
    gap_span = member.meta_fields['gap_coords'].Span()
    if (gap_seq == member.event_seq.lower() and
        dup_span == gap_span): #{
      return False
    #} end if
    gap_seq_rc = ReverseComplement(gap_seq)
    if (gap_seq_rc in ctg_seq or gap_seq_rc in t_seq): #{
      DebugMsg(self, "  Skipping microinversion gap: %s" % gap_seq_rc)
      return True
    #} end if
    return False
  #} end def

  def CheckStartAligns(self, start_aligns, gap_end, target, strand): #{
    # check whether any of the alignments looks like a local inversion type
    for align in start_aligns: #{
      if (target == align.target): #{
        DebugMsg(self, "Checking alignment: %s" % align.ToString())
        if (strand != align.strand and 0 == align.num_query_gaps and
            gap_end <= align.ctg_end): #{
          return True
        #} end if
      #} end if
    #} end for
    # none of the start alignments look like a local inversion
    return False
  #} end def

  def CheckEndAligns(self, end_aligns, gap_start, target, strand): #{
    # check whether any of the alignments looks like a local inversion type
    for align in end_aligns: #{
      if (target == align.target): #{
        DebugMsg(self, "Checking alignment: %s" % align.ToString())
        if (strand != align.strand and 0 == align.num_query_gaps and
            align.ctg_start <= gap_start): #{
          return True
        #} end if
      #} end if
    #} end for
    # none of the end alignments look like a local inversion
    return False
  #} end def

  def CheckLocalInversion(self, gap_start, gap_end, realigned_contig,
      target, strand): #{
    if (MAX_LOCAL_INVERSION > gap_start): #{
      return self.CheckStartAligns(realigned_contig.start_aligns, gap_end,
        target, strand)
    #} end if
    if (MAX_LOCAL_INVERSION > (realigned_contig.ctg_length - gap_end)): #{
      return self.CheckEndAligns(realigned_contig.end_aligns, gap_start,
        target, strand)
    #} end if
    return False
  #} end def
#} end class
