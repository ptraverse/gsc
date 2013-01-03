#! /usr/bin/env python
"""
split_candidate.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import custom modules
from utils.general import AddChr
from gene_overlap import GenesOverlapped, NearbyGenes
from alignment_functions import FixAlign, CalcAlignOverlap

# CONSTANTS
MAX_LOCAL_DIST = 5000

#Class to handle split alignments from a single contig
class SplitCandidateCls: #{
  # alignment1 in pair
  align1 = None
  # alignment2 in pair
  align2 = None
  # contig id
  contig = None
  # bp overlap in contig alignments
  ctg_overlap = None
  # fraction of shorter alignment that is overlapped by longer
  # alignment in contig coordinates
  ctg_overlap_fraction = None
  # bp's represented by the two alignments
  contig_rep = None
  # fraction of total contig represented by the two alignments
  contig_rep_fraction = None
  # alignment topology:
  #   interchromosomal, local inversion, intrachromosomal opposite-strand,
  #   junction duplication, non-colinear, read-through,
  #   intrachromosomal same-strand, end duplication
  topology = None
  # how many of the alignments in the pair could be multi-mapped
  multi_mapped = None

  def __init__(self, a1, a2, add_gene_annotation, gapped_event_found="N"): #{
    #Read in the alignment pairs (ordered by query start co-ordinate)
    a1 = FixAlign(a1)
    a2 = FixAlign(a2)
    if (a1.qstart <= a2.qstart): #{
      self.align1 = a1
      self.align2 = a2
    else: # a1.qstart > a2.qstart
      self.align1 = a2
      self.align2 = a1
    #} end if
    #attributes of the contig
    self.contig = a1.query
    self.contig_size = int(a1.query_len)
    # calculate the contig overlap and fraction
    (self.ctg_overlap, self.ctg_overlap_fraction) = (
      CalcAlignOverlap(self.align1, self.align2))
    # calculate the contig representation and fraction
    (self.contig_rep, self.contig_rep_fraction) = (
      self.CalcCtgRepresentation(self.align1, self.align2, self.ctg_overlap))
    # calculate the breakpoint distance
    self.breakpoint_distance = self.CalcBreakpointDistance()
    # calculate the target distance
    self.target_distance = self.CalcTargetDistance()
    # set the topology to unknown by default
    self.topology = "unknown"
    # count how many of the alignments are have query gaps
    self.qgaps = 0
    if (a1.qgap): #{
      self.qgaps += 1
    #} end if
    if (a2.qgap): #{
      self.qgaps += 1
    #} end if
    # count how many of the alignments are spliced
    self.spliced = 0
    if (a1.spliced): #{
      self.spliced += 1
    #} end if
    if (a2.spliced): #{
      self.spliced += 1
    #} end if
    # assume that a gapped alignment was not found
    self.gapped_event_found = gapped_event_found
    # assume that neither of the alignments can be multi-mapped
    self.multi_mapped = 0
    self.add_gene_annotation = add_gene_annotation
  #} end def

  # calculate the amount of the contig that is represented by the combination
  # of the two alignments and the fraction of the total contig length that
  # this represents
  def CalcCtgRepresentation(self, align1, align2, olap): #{
    covered_by_align1 = (align1.qend - align1.qstart) + 1
    covered_by_align2 = (align2.qend - align2.qstart) + 1
    overlap = max(0, olap)
    total_covered = covered_by_align1 + covered_by_align2 - overlap
    fraction_covered = float(total_covered) / float(align1.query_len)
    return (total_covered, fraction_covered)
  #} end def

  # label the alignment topology of the candidate
  def LabelTopology(self, min_end_dup_fract): #{
    # check whether the alignments are to the same chromosome
    if(self.align1.target == self.align2.target): #{
      # check whether the alignments are to the same strand
      if(self.align1.query_strand == self.align2.query_strand): #{
        # check whether one target region is a subset of the other
        if (self.SubsetTarget(min_end_dup_fract)): #{
          self.topology = "end-duplication"
        # check whether the query alignment order is the opposite of the
        # target alignment order
        elif("forward" != self.GetTargetOrder()): #{
          # check whether the target regions are too far apart to be
          # considered a duplication
          if (MAX_LOCAL_DIST < self.breakpoint_distance): #{
            self.topology = "intrachr-non-colinear"
          else:
            # they might be an exon duplication event
            self.topology = "junction-duplication"
          #} end if
        else:
          if (MAX_LOCAL_DIST < self.breakpoint_distance): #{
            self.topology = "intrachr-same-strand"
          else:
            self.topology = "read-through"
          #} end if
        #} end if
      # if the alignments are on opposite strands
      else:
        if (MAX_LOCAL_DIST < self.breakpoint_distance): #{
          self.topology = "intrachr-opp-strand"
        else:
          # they might be a local inversion event
          self.topology = "local-inversion"
        #} end if
      #} end if
    else:
      # they are some sort of interchromosomal chimeric event
      self.topology = "interchr"
    #} end if
  #} end def

  def CheckMultiMapping(self): #{
    if (self.align1.multi_mapped): #{
      self.multi_mapped += 1
    #} end if
    if (self.align2.multi_mapped): #{
      self.multi_mapped += 1
    #} end if
  #} end def

  def GetTargetOrder(self): #{
    if (self.align1.target != self.align2.target or
        self.align1.query_strand != self.align2.query_strand):
      return "N/A"
    #} end if
    target_order = "equal"
    # use the strand to determine the order of the alignments in
    # target (genomic) co-ordinates
    if ("+" == self.align1.query_strand): #{
      if (self.align1.tstart < self.align2.tstart): #{
        target_order = "forward"
      elif (self.align1.tstart > self.align2.tstart): #{
        target_order = "reverse"
      #} end if
    else: # "-" == self.align1.query_strand
      if (self.align1.tstart < self.align2.tstart): #{
        target_order = "reverse"
      elif (self.align1.tstart > self.align2.tstart): #{
        target_order = "forward"
      #} end if
    #} end if
    return target_order
  #} end def

  def SubsetTarget(self, min_end_dup_fract): #{
    if (self.align1.target != self.align2.target): #{
      return False
    #} end if
    start1 = min(self.align1.tstart, self.align1.tend)
    end1   = max(self.align1.tstart, self.align1.tend)
    len1 = end1 - start1 + 1
    start2 = min(self.align2.tstart, self.align2.tend)
    end2   = max(self.align2.tstart, self.align2.tend)
    len2 = end2 - start2 + 1
    if (start1 < start2): #{
      overlap = end1 - start2 + 1
    else:
      overlap = end2 - start1 + 1
    #} end if
    overlap_fraction = 0
    if (0 < overlap): #{
      min_len = min(len1, len2)
      overlap_fraction = float(overlap) / float(min_len)
    #} end if
    #if ((start1 <= start2 and end1 >= end2) or
    #    (start2 <= start1 and end2 >= end1)):
    if (min_end_dup_fract < overlap_fraction): #{
      return True
    #} end if
    return False
  #} end def

  def CalcBreakpointDistance(self): #{
    if (self.align1.target != self.align2.target): #{
      return "N/A"
    #} end if
    # determine breakpoint co-ordinate of the first block
    bp1 = self.align1.tend
    # determine breakpoint co-ordinate of the second block
    bp2 = self.align2.tstart
    # calculate the distance
    bp_dist = abs(bp1 - bp2)
    return bp_dist
  #} end def

  def CalcTargetDistance(self): #{
    if (self.align1.target != self.align2.target): #{
      return "N/A"
    #} end if
    # sort target co-ordinates of the first block
    start1 = min(self.align1.tstart, self.align1.tend)
    end1   = max(self.align1.tstart, self.align1.tend)
    # sort target co-ordinates of the second block
    start2 = min(self.align2.tstart, self.align2.tend)
    end2   = max(self.align2.tstart, self.align2.tend)
    # sort the blocks by target (genomic) start co-ordinate
    if (start1 < start2): #{
      start1_s = start1
      end1_s = end1
      start2_s = start2
      end2_s = end2
    else:
      start1_s = start2
      end1_s = end2
      start2_s = start1
      end2_s = end1
    #} end if
    td = start2_s - end1_s
    return td
  #} end def

  def BreakpointsString(self): #{
    # Determine the "breakpoint_from" target co-ordinate
    breakpoint_from = self.align1.tend
    # determine which side of the breakpoint block one represents
    if (self.align1.tstart < self.align1.tend): #{
      dir_from = "up"
    else:
      dir_from = "down"
    #} end if

    # Determine the "breakpoint_to" target co-ordinate
    breakpoint_to = self.align2.tstart
    # determine which side of the breakpoint block two represents
    if (self.align2.tstart < self.align2.tend): #{
      dir_to = "down"
    else:
      dir_to = "up"
    #} end if

    return "%s:%i(%s)-%s:%i(%s)" % (
      AddChr(self.align1.target), breakpoint_from, dir_from,
      AddChr(self.align2.target), breakpoint_to,   dir_to)
  #} end def

  def Details(self): #{
    #function formats the output lines from the alignment pair
    ctg_fields = self.contig.split(" ")
    ctg = "%s(%ibp)" % (ctg_fields[0], self.contig_size)
    target = "%s:%i-%i,%s:%i-%i" % (
      AddChr(self.align1.target), self.align1.tstart, self.align1.tend,
      AddChr(self.align2.target), self.align2.tstart, self.align2.tend)
    query = "%i-%i,%i-%i" % (
      self.align1.qstart, self.align1.qend,
      self.align2.qstart, self.align2.qend)
    align_len1 = float(self.align1.qend) - float(self.align1.qstart) + 1
    align_len2 = float(self.align2.qend) - float(self.align2.qstart) + 1
    align_frac1 = align_len1/float(self.align1.query_len)
    align_frac2 = align_len2/float(self.align2.query_len)
    align_metrics = ",".join([
      "I1:%.2f"  % float(self.align1.identity),
      "I2:%.2f"  % float(self.align2.identity),
      "AF1:%.2f" % align_frac1,
      "AF2:%.2f" % align_frac2])
    blocks1 = ",".join(["-".join(map(str, block)) for
      block in self.align1.blocks])
    blocks2 = ",".join(["-".join(map(str, block)) for
      block in self.align2.blocks])
    blocks = ";".join([blocks1, blocks2])
    olap_genes1 = GenesOverlapped(self.align1, self.add_gene_annotation)
    olap_genes2 = GenesOverlapped(self.align2, self.add_gene_annotation)
    olap_genes = ";".join([olap_genes1, olap_genes2])
    nearby_genes1 = NearbyGenes(self.align1, self.add_gene_annotation)
    nearby_genes2 = NearbyGenes(self.align2, self.add_gene_annotation)
    nearby_genes = ";".join([nearby_genes1, nearby_genes2])
    split_metrics = ",".join([
      "CO:%ibp(%.3f)" % (self.ctg_overlap, self.ctg_overlap_fraction),
      "CR:%.3f"       % self.contig_rep_fraction,
      "BD:%s"         % self.breakpoint_distance,
      "TD:%s"         % self.target_distance,
      "QG:%i"         % self.qgaps,
      "SP:%i"         % self.spliced,
      "GF:%s"         % self.gapped_event_found,
      "MM:%i"         % self.multi_mapped,
    ])
    detail_string = " ".join([
      "CTG:%s"      % ctg,
      "TOPOLOGY:%s" % self.topology,
      "TARGET:%s"   % target,
      "CONTIG:%s"   % query,
      "BREAKPOINTS:%s" % self.BreakpointsString(),
      "%s"          % align_metrics,
      "BLOCKS:%s"   % blocks,
      "GENES:%s"    % olap_genes,
      "NEARBY:%s"   % nearby_genes,
      "META:%s"     % split_metrics])
    return detail_string
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class SplitAlignmentError(Exception): #{
  def __init__(self, msg): #{
    self.msg = msg
  def __str__(self): #{
    return repr(self.msg)
#} end class
