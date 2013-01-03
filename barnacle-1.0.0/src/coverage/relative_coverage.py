#! /usr/bin/env python
"""
relative_coverage.py

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
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
    CheckConfigCommands, StrJoin, ReverseComplement, NormalizeChrID, Flatten)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, EnsureDirectoryExists,
  EnsureAbsPath, GetOutDir, CheckNewFilePath, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls
from support.samtools import SAMToolsCls
from common.coord_pair import (CoordPairCls, CutOrExtendBlocks)
from common.breakpoint import (BreakpointCls, SortBreakpoints,
  SortBreakpointTuple)

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"
SUPPORT_TYPE_KEYS = ["DW","T1","T2"]

RGB_RED  = "255,0,0"
RGB_BLUE = "0,0,255"

class ExpressionEstimatorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    CheckConfigCommands(self, "samtools")
    # create samtools object and check whether to use "chr" in chromosome IDs
    self.samtools = SAMToolsCls(self.options.p2g_path, self.options,
      log_info=self.log_info)
    self.samtools.DetermineChrUse()
    # create output file
    self.outfile = FileBoxCls(self.options.outpath, "w",
      "cannot create output file.")
    header_data =  ["ID", "C"]
    #header_data.extend(["%s%s\tC/%s%s" % (s_key, rkey, s_key, rkey) for
    #  s_key in "WT"] for rkey in ["a", "b", "min", "max", "avg"])
    #self.outfile.WriteLine(StrJoin("\t", Flatten(header_data)))
    for rkey in ["a", "b", "min", "max", "avg"]: #{
      header_data.extend("%s%s\tC/%s%s" % (s_key, rkey, s_key, rkey) for
        s_key in "WT")
    #} end for
    self.outfile.WriteLine(StrJoin("\t", header_data))
  #} end def

  def __del__(self): #{
    if (hasattr(self, "outfile") and None != self.outfile): #{
      self.outfile.Close()
    #} end if
    CloseLogFile(self)
  #} end def

  def Run(self): #{
    LogMsg(self, "Estimating relative coverage...")
    start = time.time()
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    for group in parser: #{
      DebugMsg(self, "Processing group %i" % group.id)
      # search_regions is a two-item list of CoordPairCls objects
      # with an extra "chrom" data member each
      search_regions = self.GetSearchRegions(group)
      self.CalculateCoverage(group, search_regions)
      self.WriteGroupRelativeExpression(group)
    #} end for
    LogMsg(self, "Time spent estimating relative coverage: %s" %
      TimeSpent(start))
  #} end def

  def GetSearchRegions(self, group): #{
    DebugMsg(self, "  Getting group search regions...")
    start = time.time()
    search_regions = (CoordPairCls(), CoordPairCls())
    for member in group.members: #{
      DebugMsg(self, "    Processing member %s (%s)" %
          (member.candidate_id, member.contig_info.id))
      self.GetMemberRegions(member)
      DebugMsg(self, "    Member search regions:\n%s" %
          RegionPairString(member.regions, member.IDString()))
      for (s_region, m_region) in zip(search_regions, member.regions): #{
        # get chromosomes from first member
        if (not hasattr(s_region, "chrom")): #{
          s_region.chrom = m_region.chrom
        #} end if
        s_region.Union(m_region)
      #} end for
      DebugMsg(self, "    Updated group search regions:\n%s" %
          RegionPairString(search_regions, "G%i" % group.id))
    #} end for
    DebugMsg(self, "  Time spent getting group search regions: %s" %
      TimeSpent(start))
    return search_regions
  #} end def

  def GetMemberRegions(self, member): #{
    #member.t2_support = False
    member.regions = (CoordPairCls(), CoordPairCls())
    if (member.gap): #{
      self.PreprocessGapCandidateBlocks(member)
    #} end if
    # associate contig-to-genome alignment blocks with breakpoints
    # so that they can be sorted together (to keep all members consistent)
    for (breakpoint, block_list) in zip(member.breakpoints, member.blocks): #{
      breakpoint.full_blocks = block_list
    #} end for
    # sort the breakpoints
    sorted_breakpoints = SortBreakpointTuple(member.breakpoints)
    target_length = 2 * self.options.read_length
    if ("ctg_overlap" in member.meta_fields and
        0 < member.meta_fields["ctg_overlap"]): #{
      target_length += member.meta_fields["ctg_overlap"]
    #} end if
    for (breakpoint, m_region, key) in zip(sorted_breakpoints,
        member.regions, "AB"): #{
      # cut the blocks list to the appropriate length
      ExtremeDebugMsg(self, "About to cut blocks for region %s" % key)
      m_region.target_blocks = CutOrExtendBlocks(breakpoint.full_blocks,
        target_length, from_left=("down" == breakpoint.dir.lower()))
      m_region.chrom = breakpoint.chr
      # since target blocks are sorted in forward or reverse order only need
      # to update the member search region with the first and last blocks
      m_region.Union(m_region.target_blocks[0])
      m_region.Union(m_region.target_blocks[-1])
      DebugMsg(self, "    Breakpoint: %s" % breakpoint)
      DebugMsg(self, "      Full: %s" % ",".join(("%i-%i" % block for block in
        breakpoint.full_blocks)))
      DebugMsg(self, "%s" % BlockListBEDString(m_region.chrom,
        breakpoint.full_blocks, "%s_full_blocks%s" % (member.IDString(), key),
        RGB_RED))
      #DebugMsg(self, "      Cut:  %s" % ",".join(("%i-%i" % block for block in
      #  m_region.target_blocks)))
      DebugMsg(self, "      Cut:  %s" % StrJoin(",", m_region.target_blocks))
      DebugMsg(self, "%s" % BlockListBEDString(m_region.chrom,
        m_region.target_blocks, "%s_cut_blocks%s" % (member.IDString(),
        key), RGB_BLUE))
      # create support window lists for each of the cut blocks
      DebugMsg(self, "    Creating support window lists...")
      for block in m_region.target_blocks: #{
        block.support = dict(((key, [0 for i in range(len(block))]) for
          key in SUPPORT_TYPE_KEYS))
        ExtremeDebugMsg(self, "      Block: %s (%i); %s" % (block, len(block),
          #"; ".join(("%s: %s" % (key, StrJoin(",", block.support[key])) for
          "; ".join(("%s: %i" % (key, len(block.support[key])) for
          key in SUPPORT_TYPE_KEYS))))
      #} end for
    #} end for
  #} end def

  def PreprocessGapCandidateBlocks(self, member): #{
    if (member.GapIsInternal()): #{
      # we need to cut blocks_A into two lists, by splitting the block
      # that contains the gap position (be careful of negative strand)
      DebugMsg(self, "      Before Gap: %i, After Gap: %i\n"
          "      Full blocks: %s" % (member.meta_fields['bg_end'],
          member.meta_fields['ag_start'], ",".join(("%i-%i" % block for
          block in member.blocks_A))))
      new_blocks = (list(), list())
      for block in member.blocks_A: #{
        DebugMsg(self, "      Block: %i-%i" % block)
        if (EndsBeforeGap(block, member)): #{
          # add to the blocks before the gap
          DebugMsg(self, "        Adding to before gap!")
          new_blocks[0].append(block)
        elif (StartsBeforeGap(block, member)):
          # split and add to both block lists
          DebugMsg(self, "        Splitting!")
          (before_gap, after_gap) = SplitBlock(block, member)
          DebugMsg(self, "        Before: %i-%i" % before_gap)
          DebugMsg(self, "        After:  %i-%i" % after_gap)
          new_blocks[0].append(before_gap)
          new_blocks[1].append(after_gap)
        else:
          # add to the blocks after the gap
          DebugMsg(self, "        Adding to after gap!")
          new_blocks[1].append(block)
        #} end if
      #} end for
      if (member.align_info_A.IsPosStrand()): #{
        member.blocks = tuple(reversed(new_blocks))
      else:
        member.blocks = new_blocks
      #} end if
    else:
      alignA = member.align_info_A
      alignB = member.align_info_B
      coordsA = CoordPairCls(alignA.ctg_start, alignA.ctg_end)
      coordsB = CoordPairCls(member.meta_fields['gap_start'],
          member.meta_fields['gap_end'])
      both_coords = (coordsA, coordsB)
      # fix the edge-gap (or full-contig duplication) breakpoints
      dirs = ["up", "down"]
      if (CoordsAreSorted(both_coords)): #{
        if (alignA.IsNegStrand()): #{
          dirs.reverse()
        #} end if
        breakpointA = BreakpointCls("%s:%s(%s)" % (alignA.chrom,
          alignA.genome_end, dirs[0]))
        breakpointB = BreakpointCls("%s:%s(%s)" % (alignA.chrom,
          alignB.genome_start, dirs[1]))
      else:
        if (alignA.IsPosStrand()): #{
          dirs.reverse()
        #} end if
        breakpointA = BreakpointCls("%s:%s(%s)" % (alignA.chrom,
          alignA.genome_start, dirs[0]))
        breakpointB = BreakpointCls("%s:%s(%s)" % (alignA.chrom,
          alignB.genome_end, dirs[1]))
      #} end if
      member.breakpoints = SortBreakpoints(breakpointA, breakpointB)
      if ("down" != member.breakpoints[0].dir or
          "up" != member.breakpoints[1].dir): #{
        raise ExpressionEstimatorError("incorrect breakpoint directions: "
          "%s-%s" % (member.breakpoints[0].dir, member.breakpoints[1].dir))
      # blocks are already correct, just need to make sure that they are
      # ordered properly, since gap-candidate breakpoints are always
      # left-right
      for (coords, key) in zip(both_coords, "AB"): #{
        DebugMsg(self, "      ALIGN_%s %s" % (key, coords))
      #} end for
      if ((alignA.IsPosStrand() and CoordsAreSorted(both_coords)) or
          (alignA.IsNegStrand() and not CoordsAreSorted(both_coords))): #{
        DebugMsg(self, "      Swapping blocks!")
        member.blocks = (member.blocks_B, member.blocks_A)
      #} end if
    #} end if
  #} end def

  def CalculateCoverage(self, group, search_regions): #{
    DebugMsg(self, "  Calculating group coverage...")
    start = time.time()
    # if search regions overlap, union them and only call SAMtools once
    if (search_regions[0].chrom == search_regions[1].chrom and
        search_regions[0].Overlaps(search_regions[1])): #{
      combined_region = CoordPairCls(search_regions[0])
      combined_region.Union(search_regions[1])
      combined_region.chrom = search_regions[0].chrom
      DebugMsg(self, "    Combined search regions:\n%s" %
          RegionString(combined_region, "G%i" % group.id, "C"))
      self.ProcessSearchRegion(group, combined_region, "AB")
    else:
      for (region, key) in zip(search_regions, "AB"): #{
        self.ProcessSearchRegion(group, region, key)
      #} end for
    #} end if
    self.SumSupportTypes(group)
    DebugMsg(self, "  Time spent calculating group coverage: %s" %
      TimeSpent(start))
  #} end def

  def ProcessSearchRegion(self, group, search_region, keys_to_process): #{
    DebugMsg(self, "  Processing search region(s) %s" % keys_to_process)
    count = 0
    region_str = "%s:%s-%s" % (search_region.chrom,
        search_region.min + self.options.min_overlap,
        search_region.max - self.options.min_overlap)
    for read in self.samtools.GetParsedReads(region_str,
        require_perfect=self.options.require_perfect,
        max_edit_distance = self.options.max_edit_distance): #{
      ConvertCigarToGenomeBlocks(read)
      # cut "min-overlap" off of each end of the reads alignment blocks,
      #   to ensure that support is only added to positions overlapped by
      #   at least min-overlap bases
      cut_once = CutOrExtendBlocks(read.gblocks,
          read.match_len - self.options.min_overlap, from_left=True)
      read.cut_blocks = CutOrExtendBlocks(cut_once,
          read.match_len - 2*self.options.min_overlap, from_left=False)
      #support_coords = CoordPairCls(read.pos + self.options.min_overlap,
      #  read.end_pos - self.options.min_overlap)
      support_coords = CoordPairCls(read.cut_blocks[0].min,
          read.cut_blocks[-1].max)
      if (not support_coords.Overlaps(search_region)): #{
        DebugMsg(self, "Skipping read without sufficient overlap of search "
          "region. Read: %s %i-%i (cut: %s). Search region: %s-%s" %
          (read.qname, read.pos, read.end_pos, support_coords,
          search_region.min, search_region.max))
        continue
      #} end if
      if (5 > count): #{
        if (read.reverse_strand()): #{
          read_seq = ReverseComplement(read.seq)
        else:
          read_seq = read.seq
        #} end if
        DebugMsg(self, "  Read: %s %s, CIGAR: %s %s\n    Blocks: %s\n"
          "    Cut:    %s" % (read.qname, read.coords, read.cigar,
          read_seq, StrJoin(",", read.gblocks), StrJoin(",", read.cut_blocks)))
      #} end if
      for member in group.members: #{
        self.DetermineReadSupportType(read, group.event_type,
            member.breakpoints, count)
        #if ("T2" == read.support_type): #{
        #  member.t2_support = True
        #} end if
        for (r_key, m_region) in zip("AB", member.regions): #{
          if (r_key in keys_to_process and
              support_coords.Overlaps(m_region)): #{
            self.ProcessReadWithMemberRegion(read, m_region,
                count)
          #} end if
        #} end for
      #} end for
      count += 1
      if (400 <= count): #{
        DebugMsg(self, "-"*10)
        count = 0
      #} end if
    #} end for
  #} end def

  # count is only used for debug purposes
  def DetermineReadSupportType(self, read, etype, breakpoints, count=6): #{
    if (5 > count): #{
      DebugMsg(self, "    Determining read support type: %s\n      Event "
        "type: %s, Breakpoints: %s" % (read.coords, etype, StrJoin(" / ",
        breakpoints)))
    #} end if
    if ("fusion" == etype): #{
      read.support_type = "T1"
      read.chr = NormalizeChrID(read.rname)
      for breakpoint in breakpoints: #{
        if (read.chr == NormalizeChrID(breakpoint.chr) and
            read.coords.Contains(breakpoint.coord)): #{
          read.support_type = "DW"
        #} end if
      #} end for
    elif (etype in ["ptd", "itd"]):
      breakpoint_coords = CoordPairCls([bp.coord for bp in breakpoints])
      if (read.coords.Contains(breakpoint_coords)): #{
        read.support_type = "DW"
      elif (breakpoint_coords.Contains(read.coords)):
        read.support_type = "T2"
      else:
        read.support_type = "T1"
      #} end if
    else:
      raise ExpressionEstimatorError("unrecognized event type: %s" % etype)
    #} end if
    if (5 > count): #{
      DebugMsg(self, "      Support type: %s" % (read.support_type))
    #} end if
  #} end def

  # count is only used for debug purposes
  def ProcessReadWithMemberRegion(self, read, m_region, count=6): #{
    if (5 > count): #{
      DebugMsg(self, "    Processing member region: %s, blocks: %s" %
          (m_region, StrJoin(",", m_region.target_blocks)))
    #} end if
    m_block_iter = iter(m_region.target_blocks)
    m_block = m_block_iter.next()
    for r_block in read.cut_blocks: #{
      try: #{
        while (m_block.max < r_block.min): #{
          m_block = m_block_iter.next()
        #} end while
        while (m_block.max <= r_block.max): #{
          # add support to the appropriate windows of m_block
          self.AddSupportToMemberBlock(r_block, m_block, read.support_type,
              count)
          m_block = m_block_iter.next()
        #} end while
      except StopIteration:
        break
      #} end try
      if (r_block.Overlaps(m_block)): #{
        self.AddSupportToMemberBlock(r_block, m_block, read.support_type,
            count)
      #} end if
    #} end for
    if (5 > count and self.log_info['extreme_debug']): #{
      for s_key in SUPPORT_TYPE_KEYS: #{
        DebugMsg(self, "      %s: %s" % (s_key, StrJoin("; ", (StrJoin(",",
          ("%2i" % s for s in block.support[s_key])) for block in
          m_region.target_blocks))))
      #} end for
    #} end if
  #} end def

  def AddSupportToMemberBlock(self, r_block, m_block, support_type, count=6): #{
    supported_indices = m_block.Copy()
    supported_indices.Intersect(r_block)
    if (5 > count): #{
      DebugMsg(self, "      Adding support from read block %s, to member "
        "block %s. Intersection: %s" % (r_block, m_block, supported_indices))
    #} end if
    supported_indices.MoveMin(-m_block.min)
    supported_indices.MoveMax(-m_block.min)
    if (5 > count): #{
      DebugMsg(self, "      Supported windows: %s" % supported_indices)
    #} end if
    for i in range(supported_indices.min, supported_indices.max+1): #{
      m_block.support[support_type][i] += 1
    #} end for
  #} end def

  def SumSupportTypes(self, group): #{
    for member in group.members: #{
      DebugMsg(self, "  MEMBER: %s %s" % (member.IDString(),
        member.read_to_ctg_unique))
      for (m_region, key) in zip(member.regions, "AB"): #{
        DebugMsg(self, "    REGION %s: chr%s %s" % (key, m_region.chrom,
          m_region))
        m_region.max_vals = {"T": set(), "W": -2*member.read_to_ctg_unique}
        for block in m_region.target_blocks: #{
          block.support["T"] = list()
          block.support["W"] = list()
          block.max_vals = {"T": set(), "W": -2*member.read_to_ctg_unique}
          #support_lists = zip((block.support[s_key] for s_key in
          #    SUPPORT_TYPE_KEYS))
          #for (dw_sup, t1_sup, t2_sup) in support_lists: #{}
          num_wins = min(len(block.support[s_key]) for s_key in
            SUPPORT_TYPE_KEYS)
          for i in xrange(num_wins): #{
            (dw_sup, t1_sup, t2_sup) = (block.support[s_key][i] for s_key in
                SUPPORT_TYPE_KEYS)
            t_sup = dw_sup + t1_sup + t2_sup
            block.support["T"].append(t_sup)
            w_sup = t_sup - member.read_to_ctg_unique
            if (0 < t2_sup): #{
              w_sup -= 2 * member.read_to_ctg_unique
            #} end if
            block.support["W"].append(w_sup)
            if (w_sup > block.max_vals["W"]): #{
              block.max_vals["W"] = w_sup
              #ExtremeDebugMsg(self, "Resetting block.max_vals[\"T\"]: "
              #  "w_sup = %i, max_w = %i" % (w_sup, block.max_vals["W"]))
              block.max_vals["T"] = set()
            #} end if
            if (w_sup == block.max_vals["W"]): #{
              #ExtremeDebugMsg(self, "Updating block.max_vals[\"T\"]: "
              #  "w_sup = %i, max_w = %i" % (w_sup, block.max_vals["W"]))
              block.max_vals["T"].add(t_sup)
            #} end if
          #} end for
          if (block.max_vals["W"] > m_region.max_vals["W"]): #{
            m_region.max_vals["W"] = block.max_vals["W"]
            #ExtremeDebugMsg(self, "Resetting m_region.max_vals[\"T\"]: "
            #  "block_max = %i, m_region_max = %i" % (block.max_vals["W"],
            #  m_region.max_vals["W"]))
            m_region.max_vals["T"] = set()
          #} end if
          if (block.max_vals["W"] == m_region.max_vals["W"]): #{
            #ExtremeDebugMsg(self, "Updating m_region.max_vals[\"T\"]: "
            #  "block_max = %i, m_region_max = %i" % (block.max_vals["W"],
            #  m_region.max_vals["W"]))
            m_region.max_vals["T"].update(block.max_vals["T"])
          #} end if
        #} end for
        for s_key in (SUPPORT_TYPE_KEYS + ["T", "W"]): #{
          DebugMsg(self, "      %2s: %s" % (s_key, StrJoin("; ", (StrJoin(",",
            ("%3i" % s for s in block.support[s_key])) for block in
            m_region.target_blocks))))
        #} end for
        m_region.max_vals["T"]= max(m_region.max_vals["T"])
        DebugMsg(self, "      Max W: %i, T: %i" % (m_region.max_vals["W"],
          m_region.max_vals["T"]))
      #} end for
    #} end for
  #} end def

  def WriteGroupRelativeExpression(self, group): #{
    min_member = None
    max_member = None
    sums = {'C':0, 'Wa':0, 'Ta':0, 'Wb':0, 'Tb':0}
    for member in group.members: #{
      if (None == min_member or
          min_member.read_to_ctg_unique > member.read_to_ctg_unique): #{
        min_member = member
      #} end if
      if (None == max_member or
          max_member.read_to_ctg_unique < member.read_to_ctg_unique): #{
        max_member = member
      #} end if
      sums["C"] += member.read_to_ctg_unique
      for (m_region, rkey) in zip(member.regions, "ab"): #{
        for s_key in "WT": #{
          sums["%s%s" % (s_key,rkey)] += m_region.max_vals[s_key]
        #} end for
      #} end for
      self.WriteMemberRelativeExpression(member.IDString(),
        member.read_to_ctg_unique, ([m_region.max_vals[s_key] for
        s_key in "WT"] for m_region in member.regions))
        #member.regions[0].max_vals["W"], member.regions[0].max_vals["T"],
        #member.regions[1].max_vals["W"], member.regions[1].max_vals["T"])
    #} end for
    self.WriteMemberRelativeExpression("%imin" % group.id,
      min_member.read_to_ctg_unique, ([m_region.max_vals[s_key] for
        s_key in "WT"] for m_region in min_member.regions))
      #min_member.regions[0].max_vals["W"], min_member.regions[0].max_vals["T"],
      #min_member.regions[1].max_vals["W"], min_member.regions[1].max_vals["T"])
    self.WriteMemberRelativeExpression("%imax" % group.id,
      max_member.read_to_ctg_unique, ([m_region.max_vals[s_key] for
        s_key in "WT"] for m_region in max_member.regions))
    avg_vals = dict((s_key, cov/float(len(group.members))) for (s_key,cov) in
        sums.iteritems())
    self.WriteMemberRelativeExpression("%iavg" % group.id,
      avg_vals["C"], ([avg_vals["%s%s" % (s_key,rkey)] for s_key in "WT"] for
        rkey in "ab"))
  #} end def

  def WriteMemberRelativeExpression(self, id, c_cov, w_and_t_coverage): #{
    data_list = list()
    data_list.append(id)
    data_list.append(c_cov)
    (a_cov_pair, b_cov_pair) = w_and_t_coverage
    data_list.extend(CoveragePairOutput(c_cov, a_cov_pair))
    data_list.extend(CoveragePairOutput(c_cov, b_cov_pair))
    min_cov_pair = map(min, zip(a_cov_pair, b_cov_pair))
    #ExtremeDebugMsg(self, "A: %s, B: %s, MIN: %s" %
    #    (a_cov_pair, b_cov_pair, min_cov_pair))
    data_list.extend(CoveragePairOutput(c_cov, min_cov_pair))
    max_cov_pair = map(max, zip(a_cov_pair, b_cov_pair))
    data_list.extend(CoveragePairOutput(c_cov, max_cov_pair))
    avg_cov_pair = ((a+b)/2. for (a,b) in zip(a_cov_pair, b_cov_pair))
    data_list.extend(CoveragePairOutput(c_cov, avg_cov_pair))
    self.outfile.WriteLine(StrJoin("\t", data_list))
  #} end def
#} end class

def BlockListBEDString(chr, raw_blocks, name, colour): #{
  block_list = [CoordPairCls(block) for block in raw_blocks]
  bed_data = list()
  # 1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or
  # scaffold (e.g. scaffold10671).
  if (not chr.startswith("chr")): #{
    chr = "chr%s" % chr
  #} end if
  bed_data.append(chr)
  # 2. chromStart - The starting position of the feature in the chromosome or
  #    scaffold. The first base in a chromosome is numbered 0.
  start = min(block_list[0].start, block_list[-1].end) - 1
  bed_data.append(start)
  # 3. chromEnd - The ending position of the feature in the chromosome or
  #    scaffold. The chromEnd base is not included in the display of the
  #    feature. For example, the first 100 bases of a chromosome are defined
  #    as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
  end = max(block_list[0].start, block_list[-1].end)
  bed_data.append(end)
  # The 9 additional optional BED fields are:
  # 4. name - Defines the name of the BED line. This label is displayed to the
  #    left of the BED line in the Genome Browser window when the track is
  #    open to full display mode or directly to the left of the item in pack
  #    mode.
  bed_data.append(name)
  # 5. score - A score between 0 and 1000. If the track line useScore
  #    attribute is set to 1 for this annotation data set, the score value
  #    will determine the level of gray in which this feature is displayed
  #    (higher numbers = darker gray).
  bed_data.append(1000)
  # 6. strand - Defines the strand - either '+' or '-'.
  if (block_list[0].start < block_list[-1].end): #{
    is_positive = True
    bed_data.append("+")
  else:
    is_positive = False
    bed_data.append("-")
  #} end if
  # 7. thickStart - The starting position at which the feature is drawn
  #    thickly (for example, the start codon in gene displays).
  bed_data.append(start)
  # 8. thickEnd - The ending position at which the feature is drawn thickly
  #    (for example, the stop codon in gene displays).
  bed_data.append(end)
  # 9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track
  #    line itemRgb attribute is set to "On", this RBG value will determine
  #    the display color of the data contained in this BED line. NOTE: It is
  #    recommended that a simple color scheme (eight colors or less) be used
  #    with this attribute to avoid overwhelming the color resources of the
  #    Genome Browser and your Internet browser.
  bed_data.append(colour)
  # 10. blockCount - The number of blocks (exons) in the BED line.
  bed_data.append(len(block_list))
  # 11. blockSizes - A comma-separated list of the block sizes. The number of
  #    items in this list should correspond to blockCount.
  block_sizes = StrJoin(",", (len(block) for block in block_list))
  bed_data.append(block_sizes)
  # 12. blockStarts - A comma-separated list of block starts. All of the
  #    blockStart positions should be calculated relative to chromStart. The
  #    number of items in this list should correspond to blockCount.
  if (is_positive): #{
    block_starts = StrJoin(",", ((block.start-1)-start for
      block in block_list))
  else:
    block_starts = StrJoin(",", ((block.end-1)-start for
      block in block_list))
  #} end if
  bed_data.append(block_starts)
  return StrJoin(" ", bed_data)
#} end def

def RegionString(region, id, key): #{
  return ("chr%s %s %s%s" % (region.chrom, region, id, key)).replace("-", " ")
#} end def

def RegionPairString(region_pair, id, indent=""): #{
  delim = "\n%s" % indent
  return (delim.join((RegionString(region, id, key) for (region, key) in
    zip(region_pair, "AB"))))
#} end def

def CoordsAreSorted(both_coords): #{
  return both_coords[0].min < both_coords[1].min
#} end def

def EndsBeforeGap(block, member): #{
  gap_coord = member.meta_fields['bg_end']
  if (member.align_info_A.IsPosStrand()): #{
    return (block[1] <= gap_coord)
  else:
    return (block[1] >= gap_coord)
  #} end if
#} end def

def StartsBeforeGap(block, member): #{
  gap_coord = member.meta_fields['bg_end']
  if (member.align_info_A.IsPosStrand()): #{
    return (block[0] <= gap_coord)
  else:
    return (block[0] >= gap_coord)
  #} end if
#} end def

def SplitBlock(block, member): #{
  gap_coord = member.meta_fields['bg_end']
  before_gap = (block[0], gap_coord)
  if (member.align_info_A.IsPosStrand()): #{
    after_gap = (gap_coord+1, block[1])
  else:
    after_gap = (gap_coord-1, block[1])
  #} end if
  return (before_gap, after_gap)
#} end def

def ConvertCigarToGenomeBlocks(read): #{
  read.gblocks = list()
  read.match_len = 0
  start = read.pos
  for cigar_field in read.cigar_list(): #{
    count = int(cigar_field[:-1])
    type = cigar_field[-1]
    # if the cigar field is an alignment match, create a block
    if ("M" == type): #{
      end = start + count - 1
      new_block = CoordPairCls(start, end)
      read.gblocks.append(new_block)
      read.match_len += len(new_block)
      start = end + 1
    # if the cigar field is a deletion from the reference or
    #   a skipped region from the reference (intron), move the start
    elif (type in "DN"): #{
      start += count
    # otherwise, the cigar field should be an insertion to the reference (I),
    #   a hard or soft clipped region (H, S), or a padding block
    elif (type not in "IHPS"): #{
      raise ExpressionEstimatorError("Unrecognized CIGAR field: %s" % type)
    #} end if
  #} end for
  read.end_pos = read.gblocks[-1].max
  read.coords = CoordPairCls(read.pos, read.end_pos)
#} end def

def CoveragePairOutput(c_cov, cov_pair): #{
  #cov_pair = tuple(cov_pair)
  #print "COV_PAIR:", cov_pair
  #print "COV_PAIR:", tuple(cov_pair)
  #print "COV_PAIR:", list(cov_pair)
  #(w_cov, t_cov) = tuple(cov_pair)
  (w_cov, t_cov) = cov_pair
  #print "GOOD"
  return CoverageOutput(c_cov, w_cov) + CoverageOutput(c_cov, t_cov)
#} end def

def CoverageOutput(c_cov, other_cov): #{
  if (0 >= other_cov): #{
    ratio = 999999
  else:
    ratio = float(c_cov)/float(other_cov)
  #} end if
  return [other_cov, ratio]
#} end def

#### EXCEPTION CLASSES ####
class ExpressionEstimatorError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Calculates maximum total coverage of exons involved "
      "in events, then compares that to event coverage to estimate relative "
      "coverage of event transcripts.")
  args = [ "LIB", "BARNACLE_FILE",  "P2G_FILE" ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-r", "--read-length",
    type="int",
    help="The length of the reads in the paired-read to genome file. "
      "[default=%default].")
  parser.add_option("--min-overlap",
    type="int", metavar="N",
    help="Require that reads overlap positions by at least N bp. "
      "[default: %default]")
  parser.add_option("--allow-mismatches",
    action="store_false", dest="require_perfect",
    help="Allow gaps and mismatches in read-to-genome alignments. [default]")
  parser.add_option("--require-perfect",
    action="store_true", dest="require_perfect",
    help="Only count reads with perfect read-to-genome alignments.")
  parser.add_option("-e", "--max-edit-distance",
    type="int", metavar="N",
    help="Only count reads with an edit distance not greater than N. "
      "[default: %default]")
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  misc_group.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  misc_group.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(read_length=75,
                      min_overlap=5,
                      max_edit_distance=1,
                      dpt=False,
                      force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (not opts_good): #{
    ErrMsg("bad option") #TODO
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  CheckFilePath(options.p2g_path, "read-to-genome alignments", path_errors)
  # get and check the output path
  options.output_dir = GetOutDir(os.path.dirname(options.barnacle_path),
    "relative_coverage")
  #CheckDirPath(options.output_dir, "output", path_errors,
  #  create=True, replace=True)
  EnsureDirectoryExists(options.output_dir)
  if (opts_good and 0 == len(path_errors)): #{
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "relative_coverage", options.output_dir)
    barnacle_file_name = os.path.basename(options.barnacle_path)
    options.outpath = os.path.join(options.output_dir,
        "%s.relative_coverage.txt" % barnacle_file_name)
    if (not options.force): #{
      CheckNewFilePath(options.outpath, "relative coverage output file",
        path_errors)
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
    options.p2g_path      = EnsureAbsPath(args[2])
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = ExpressionEstimatorCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.Run()
      except (MyError), e:
        ErrMsg("ERROR while estimating relative coverage:\n  %s" % e) #TODO
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
