#! /usr/bin/env python
"""
ptd_predictor.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
from utils.messages import LogMsg, DebugMsg
from base_predictor import BasePredictorCls

# constants
#MIN_IMPROVEMENT = 10
MIN_IMPROVE_FRACTION = 1.5

class PTDPredictorCls(BasePredictorCls): #{
  def __init__(self, options, log_info=None): #{
    BasePredictorCls.__init__(self, 'ptd', 'ptd', 'partial tandem duplication',
      options, log_info=log_info)
    self.good_topologies = set([
      'end-duplication',
      'intrachr-non-colinear',
      'junction-duplication',
    ])
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
    # if the event does not match two exon boundaries, it is not a PTD
    if (2 > member.meta_fields['exon_bounds']): #{
      DebugMsg(self, "  %s: Too few exon boundaries" % member.IDString())
      return False
    #} end if
    # if either gene set is empty, then the member is not
    # a partial tandem duplication
    if (not member.OverlapsGene("A", exclude_UTRs=True) or
        not member.OverlapsGene("B", exclude_UTRs=True)):
      DebugMsg(self, "  %s: Empty gene set!" % member.IDString())
      return False
    #} end if
    # if the gene sets overlap at all, then the event might be
    # a full exon duplication
    if (member.GeneSetsOverlap()): #{
      DebugMsg(self, "  %s: Overlapping gene sets!" % member.IDString())
      return True
    #} end if
    DebugMsg(self, "  %s: Disjoint gene sets!" % member.IDString())
    return False
  #} end def

  #def ReprocessMember(self, member, contig_alignment, contig_seq=None): #{
  def ReprocessMember(self, member, realigned_contig): #{
    contig_alignment = realigned_contig.best_align
    if (None == contig_alignment): #{
      DebugMsg(self, "  Contig %s could not be aligned to any wild-type "
          "transcript." % realigned_contig.id)
      return True
    #} end if
    if (contig_alignment.perfect): #{
      DebugMsg(self, "  False positive: perfect alignment (%i)" %
        contig_alignment.match)
      return False
    #} end if
    #new_span = (contig_alignment.ctg_end - contig_alignment.ctg_start) + 1
    #DebugMsg(self, "  New span: %i" % new_span)
    new_match = contig_alignment.match
    min_match = int(max(member.align_info_A.ContigSpan(),
      member.align_info_B.ContigSpan()) * MIN_IMPROVE_FRACTION)
    DebugMsg(self, "  New match: %i, Min: %i" % (new_match, min_match))
    # if the new alignment is significantly longer (in the contig) than
    # the split contig-to-genome alignments, then the contig is probably
    # not really a partial tandem duplication
    #if (MIN_IMPROVEMENT < (new_span - member.align_info_A.ContigSpan()) and
    #    MIN_IMPROVEMENT < (new_span - member.align_info_B.ContigSpan())): #{}
    if (min_match < new_match): #{
      DebugMsg(self, "  False positive: better alignment")
      return False
    #} end if
    DebugMsg(self, "  Good ptd prediction")
    return True
  #} end def
#} end class
