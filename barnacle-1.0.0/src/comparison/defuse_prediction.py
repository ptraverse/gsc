#! /usr/bin/env python
"""
defuse_prediction.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
from utils.general import NormalizeChrID
from utils.messages import DebugMsg
from common.breakpoint import BreakpointCls, SortBreakpoints

# CONSTANTS
DEFUSE_COMPARISON_HEADER="\t".join(["cluster_id",
  "chrom1","gene_id1","gene_name1","gene_part1","breakpoint1",
  "chrom2","gene_id2","gene_name2","gene_part2","breakpoint2",
  "probability","coverage","interchr","inversion","readthrough","deletion",
  "lib_name","sequence"])
ALIGN_RESULTS_HEADER="\t".join(["cluster_id","partners","coverage","result",
  "contig(s)"])
MATCH_THRESHOLD=5

class DefusePredictionCls: #{
  def __init__(self, data_str, field_cols, log_info=None): #{
    self.log_info = log_info
    input_dicts = self.GetInputDicts(data_str, field_cols)
    self.SetupDataFields(input_dicts[0])
    if (0 < len(input_dicts[1][0]) and 0 < len(input_dicts[1][1])): #{
      self.SetupBreakpointFields(input_dicts[1])
    #} end if
    self.quality = 0
    self.InitializeAlignResults()
  #} end def

  def GetInputDicts(self, data_str, field_cols): #{
    data_list = data_str.split()
    data_fields = dict()
    breakpoint_fields = [dict(),dict()]
    for field_name in field_cols: #{
      column = field_cols[field_name]
      if (len(data_list) <= column): #{
        raise DefusePredictionError("truncated record:\n%s\nMissing %s:%i" %
          (data_str, field_name, column))
      #} end if
      if (field_name[-1] in ["1", "2"]): #{
        bp_field_name = field_name[:-1]
        bp_index = int(field_name[-1])-1
        #DebugMsg(self, "Field: %s, Breakpoint field: %s, Breakpoint index: %i" %
        #  (field_name, bp_field_name, bp_index))
        breakpoint_fields[bp_index][bp_field_name] = data_list[column]
      else:
        data_fields[field_name] = data_list[column]
      #} end if
    #} end for
    return (data_fields, breakpoint_fields)
  #} end def

  def SetupDataFields(self, data_fields): #{
    self.id  = data_fields.get('cluster_id')
    self.seq = data_fields.get('splitr_sequence')
    if (None != self.seq): #{
      self.breakpoint = self.seq.find("|")
    else:
      self.breakpoint = None
    #} end if
    if ('splitr_count' in data_fields):
      self.coverage = int(data_fields['splitr_count'])
    else:
      self.coverage = None
    #} end if
    #readthrough_fields = [data_fields['adjacent'], data_fields['altsplice'],
    if ('read_through' in data_fields): #{
      readthrough_fields = [data_fields['read_through']]
      if ('altsplice' in data_fields): #{
        readthrough_fields.append(data_fields['altsplice'])
      #} end if
      self.readthrough = ("Y" in readthrough_fields)
    else:
      self.readthrough = None
    #} end if
    if ('deletion' in data_fields): #{
      self.deletion = ("Y" == data_fields['deletion'])
    else:
      self.deletion = None
    #} end if
    if ('interchromosomal' in data_fields): #{
      self.interchr = ("Y" == data_fields['interchromosomal'])
    else:
      self.interchr = None
    # end if
    if ('inversion' in data_fields): #{
      self.inversion = ("Y" == data_fields['inversion'])
    else:
      self.inversion = None
    #} end if
    self.lib = data_fields.get('library_name')
    if ('probability' in data_fields):
      self.prob = float(data_fields.get('probability'))
    else:
      self.prob = None
    # end if
  #} end def

  def SetupBreakpointFields(self, breakpoint_fields): #{
    breakpoint0 = DefuseBreakpointCls(breakpoint_fields[0])
    breakpoint1 = DefuseBreakpointCls(breakpoint_fields[1])
    sorted_breakpoints = SortBreakpoints(breakpoint0, breakpoint1)
    (self.breakpointA, self.breakpointB) = sorted_breakpoints
  #} end def

  def InitializeAlignResults(self): #{
    self.full_align = False
    self.part_align = False
    #self.best_full_match = 0
    #self.best_tfull_match = 0
    #self.best_part_match = 0
    self.best_aligns = {
      'full': AlignCollectionCls(),
      'tfull': AlignCollectionCls(),
      'part': AlignCollectionCls(),
    }
    self.realigned = False
    self.any_full_realign = False
    self.any_part_realign = False
    # self.full_realign[kvalue] = bool
    self.full_realign = dict()
    # self.part_realign[kvalue] = bool
    self.part_realign = dict()
    # self.best_realigns[kvalue][align_type] = AlignCollectionCls()
    self.best_realigns = dict()
  #} end def

  #def UpdateBestAligns(self, new_align): #{
  #  if (MATCH_THRESHOLD < (self.best_match - new_align.match)): #{
  #    return
  #  #} end if
  #  self.best_aligns.append(new_align)
  #  if (new_align.match > self.best_match): #{
  #    self.best_match = new_align.match
  #    new_best_aligns = list()
  #    for align in self.best_aligns: #{
  #      if (MATCH_THRESHOLD > (self.best_match - align.match)): #{
  #        new_best_aligns.append(align)
  #      #} end if
  #    #} end for
  #    self.best_aligns = new_best_aligns
  #  #} end if
  #} end def

  def ToString(self): #{
    data_list = [self.id,self.breakpointA.FullString(),self.breakpointB.FullString(),
      self.prob,self.coverage,self.interchr,self.inversion,self.readthrough,
      self.deletion,self.lib,self.seq]
    data_str = "\t".join(map(str, data_list))
    return data_str.replace("True","1").replace("False","0")
  #} end def

  def ShortString(self): #{
    return "%s: %s/%s" % (self.id, self.breakpointA.gene_name,
      self.breakpointB.gene_name)
  #} end def

  def AlignString(self): #{
    result = "no_align"
    contig_list = list()
    if (self.realigned): #{
      if (self.any_full_realign): #{
        result = "full_realign"
        contig_list.extend(self.GetRealignContigs("full"))
        contig_list.extend(self.GetRealignContigs("tfull"))
      elif (self.any_part_realign):
        result = "part_realign"
        contig_list.extend(self.GetRealignContigs("part"))
      #} end if
    else:
      if (self.full_align): #{
        result = "full_align"
        contig_list.extend(self.GetAlignContigs("full"))
        contig_list.extend(self.GetAlignContigs("tfull"))
      elif (self.part_align):
        result = "part_align"
        contig_list.extend(self.GetAlignContigs("part"))
      #} end if
    #} end if
    #contigs = ",".join(sorted((align.target for align in self.best_aligns)))
    contigs = ",".join(sorted(contig_list))
    data_list = [self.id, "%s/%s" % (self.breakpointA.gene_name,
      self.breakpointB.gene_name), str(self.coverage), result, contigs]
    data_str = "\t".join(data_list)
    return data_str
  #} end def

  def GetRealignContigs(self, align_type): #{
    contig_list = list()
    #for kvalue in self.best_realigns[align_type].keys(): #{}
    for kvalue in sorted(self.best_realigns.keys()): #{
      contig_list.extend(("k%s:%s" % (kvalue, align.target) for align in
        self.best_realigns[kvalue][align_type]))
    #} end for
    return contig_list
  #} end def

  def GetAlignContigs(self, align_type): #{
    return self.best_aligns[align_type].TargetList()
  #} end def

  def ShortSequence(self, dist_from_break): #{
    break_pos = self.breakpoint
    if (-1 == break_pos): #{
      raise DefusePredictionError("breakpoint not marked with \"|\" in "
        "sequence for Defuse prediction %i" % self.id)
    #} end if
    before_break = self.seq[break_pos-dist_from_break:break_pos]
    after_break = self.seq[break_pos+1:break_pos+dist_from_break+1]
    return before_break + after_break
  #} end def
#} end class

class AlignCollectionCls: #{
  def __init__(self): #{
    self.Reset()
  #} end def

  def Update(self, new_align): #{
    if (MATCH_THRESHOLD < (self.best_match - new_align.match)): #{
      return
    #} end if
    self.aligns.append(new_align)
    if (new_align.match > self.best_match): #{
      self.best_match = new_align.match
      new_align_list = list()
      for align in self.aligns: #{
        if (MATCH_THRESHOLD > (self.best_match - align.match)): #{
          new_align_list.append(align)
        #} end if
      #} end for
      self.best_aligns = new_align_list
    #} end if
  #} end def

  def Reset(self): #{
    self.best_match = 0
    self.aligns = list()
  #} end def

  def TargetList(self): #{
    return (align.target for align in self.aligns)
  #} end def
#} end class

class DefuseBreakpointCls(BreakpointCls): #{
  def __init__(self, data_fields): #{
    #chrom = data_fields['gene_chromosome'].upper()
    #chrom = chrom.replace("chr","").replace("MT","M")
    chrom = NormalizeChrID(data_fields['gene_chromosome'])
    breakpoint_str = "%s:%s(up)" % (chrom, data_fields['genomic_break_pos'])
    BreakpointCls.__init__(self, breakpoint_str)
    self.gene_id   = data_fields['gene']
    self.gene_part = data_fields['gene_location']
    self.gene_name = data_fields['gene_name']
  #} end def

  def FullString(self): #{
    data_list = [self.chr,self.gene_id,self.gene_name,self.gene_part,
      self.coord,]
    return "\t".join(map(str,data_list))
  #} end def

  def ToString(self): #{
    return BreakpointCls.ToString(self)
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class DefusePredictionError(MyError): #{
  pass
#} end class
