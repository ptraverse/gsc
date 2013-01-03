#! /usr/bin/env python
"""
group.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import time

# import custom modules
from utils.log import CloseLogFile
from utils.error import MyError
from utils.general import SetupMainClass, TimeSpent, AddChr, NormalizeChrID
from utils.messages import LogMsg, ExtremeDebugMsg
from utils.files_paths import FileBoxCls
from support.samtools import SAMToolsCls, SAMToolsError

# CONSTANTS
UPSTREAM_INDEX   = 0
DOWNSTREAM_INDEX = 1
SAM_BOTH         = 0
SAM_UPSTREAM     = 1
SAM_DOWNSTREAM   = 2

class P2GGroupCls: #{
  def __init__(self, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    self.regions  = list([None, None])
    self.use_chr  = False
    if (hasattr(options, "use_chr")): #{
      self.use_chr = options.use_chr
    #} end if
    #ExtremeDebugMsg(self, "Should I use chr? %s" % self.use_chr)
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def ParseGroupLine(self, group_line): #{
    (group_id, coordsA, coordsB, self.topology, self.ctg_overlap,
     self.ctg_ids) = group_line.split("\t")
    self.group_id = int(group_id)
    self.coords = list([PairToGenomeCoordsCls(coordsA),
      PairToGenomeCoordsCls(coordsB)])
    if ("N/A" != self.ctg_overlap): #{
      self.ctg_overlap = int(self.ctg_overlap)
    #} end if
  #} end def

  def GetPairToGenomeSupport(self): #{
    start_time = time.time()
    # get the length of the genomic region to search
    self.GetRegionLength()
    # get the coordinates of the upstream region
    self.GetRegionCoords(UPSTREAM_INDEX)
    # get the coordinates of the downstream region
    self.GetRegionCoords(DOWNSTREAM_INDEX)
    # determine which region to use samtools for
    self.ChooseSAMRegion()
    ExtremeDebugMsg(self, "Time spent preparing for samtools: %s" %
      TimeSpent(start_time))
    # run samtools on the chosen region(s)
    if (SAM_BOTH == self.what_to_sam): #{
      sam_path = self.RunSAMTools(self.primary_index,
          out_file="sam_out_tmp_1", sort=True)
      sam_path2 = self.RunSAMTools(self.secondary_index,
          out_file="sam_out_tmp_2", sort=True)
      self.PrepareSecondarySamResults(sam_path2)
    else:
      sam_path = self.RunSAMTools(self.primary_index)
      self.secondary_sam_file = None
    #} end if
    # parse the samtools results
    self.ParseSamToolsResults(sam_path)
  #} end def

  def GetRegionLength(self): #{
    self.region_len = 1.25 * float(self.options.frag_len)
    if ("N/A" != self.ctg_overlap and 0 < self.ctg_overlap): #{
      self.region_len += self.ctg_overlap
    #} end if
    ExtremeDebugMsg(self,
      "Contig Overlap: %s\n" % str(self.ctg_overlap) +
      "Region Length: %i" % self.region_len)
  #} end def

  def GetRegionCoords(self, region_index): #{
    region = SearchRegionCls(self.options, self.log_info)
    region.GetCoords(self.topology, self.coords[region_index], self.region_len)
    ExtremeDebugMsg(self,
      "Full Coords: %s\n" % self.coords[region_index].ToString() +
      "Region %i: %s" % (region_index, region.ToString()))
    if (region.left > region.right): #{
      LogMsg(self, "WARNING: bad region %i: \"%i-%i\". " %
        (region_index, region.left, region.right) +
        "For group %i with contigs: %s" % (self.group_id, self.ctg_ids))
    #} end if
    self.regions[region_index] = region
  #} end def

  def AdjustGapRegions(self): #{
    # downstream left and right are both the right side of the gap
    self.regions[UPSTREAM_INDEX].right  = self.regions[DOWNSTREAM_INDEX].right
    # upstream left and right are both the left side of the gap
    self.regions[DOWNSTREAM_INDEX].left = self.regions[UPSTREAM_INDEX].left
    # add the blocks information from the contig to genome alignment
    if ("duplication" in self.topology): #{
      self.regions[DOWNSTREAM_INDEX].blocks = \
        self.coords[UPSTREAM_INDEX].blocks
    else:
      self.regions[DOWNSTREAM_INDEX].blocks = self.MergeBlocks()
    #} end if
    ExtremeDebugMsg(self, "Adjusting gap regions...\n"
      "Region %i: %s\n" %
        (UPSTREAM_INDEX, self.regions[UPSTREAM_INDEX].ToString()) +
      "Region %i: %s" %
        (DOWNSTREAM_INDEX, self.regions[DOWNSTREAM_INDEX].ToString()))
  #} end def

  def MergeBlocks(self): #{
    ExtremeDebugMsg(self, "MERGING: %s; %s" %
        (self.coords[UPSTREAM_INDEX].blocks,
         self.coords[DOWNSTREAM_INDEX].blocks))
    main_blocks = self.coords[UPSTREAM_INDEX].blocks.split(",")
    gap_blocks  = self.coords[DOWNSTREAM_INDEX].blocks.split(",")
    merged_blocks = list()
    (main_index, gap_index) = (0, 0)
    while (main_index < len(main_blocks) and gap_index < len(gap_blocks)): #{
      main_left = main_blocks[main_index].split("-")[0]
      gap_left  = gap_blocks[gap_index].split("-")[0]
      if (main_left < gap_left): #{
        merged_blocks.append(main_blocks[main_index])
        main_index += 1
      else: # gap_left <= main_left
        merged_blocks.append(gap_blocks[gap_index])
        gap_index += 1
      #} end if
    #} end while
    while (main_index < len(main_blocks)): #{
      merged_blocks.append(main_blocks[main_index])
      main_index += 1
    #} end while
    while (gap_index < len(gap_blocks)): #{
      merged_blocks.append(gap_blocks[gap_index])
      gap_index += 1
    #} end while
    merged_blocks_str = ",".join(merged_blocks)
    ExtremeDebugMsg(self, "Merged blocks: %s" % merged_blocks_str)
    return merged_blocks_str
  #} end def

  def ChooseSAMRegion(self): #{
    self.what_to_sam     = SAM_BOTH
    self.primary_index   = UPSTREAM_INDEX
    self.secondary_index = DOWNSTREAM_INDEX
    # if this is a gap candidate
    if ("gap" in self.topology): #{
      # adjust the region coords
      self.AdjustGapRegions()
    # if the upstream breakpoint is on the left side
    elif ("left" == self.coords[UPSTREAM_INDEX].breakside): #{
      # sam the downstream region
      self.what_to_sam     = SAM_DOWNSTREAM
      self.primary_index   = DOWNSTREAM_INDEX
      self.secondary_index = UPSTREAM_INDEX
    # if the downstream breakpoint is on the left side
    elif ("left" == self.coords[DOWNSTREAM_INDEX].breakside): #{
      # sam the upstream region
      self.what_to_sam = SAM_UPSTREAM
    #} end if
    ExtremeDebugMsg(self, "What to SAM: %s" %
      ("both", "up", "down")[self.what_to_sam])
  #} end def

  def RunSAMTools(self, region_index, out_file=None, sort=False): #{
    start = time.time()
    ExtremeDebugMsg(self, "Running samtools...")
    samtools = SAMToolsCls(self.options.p2g_path, self.options,
      log_info=self.log_info, out_file=out_file, sort=sort, paired=True)
    try:
      samtools.Run(self.regions[region_index].ToString(self.use_chr))
    except SAMToolsError, e:
      raise P2GGroupError \
        ("Error running samtools on group %i with contig(s): %s\n%s" %
         (self.group_id, self.ctg_ids, e))
    # end try
    ExtremeDebugMsg(self, "Time spent running samtools: %s" % TimeSpent(start))
    return samtools.out_file_path
  #} end def

  def PrepareSecondarySamResults(self, sam_path): #{
    fail_msg = "cannot open secondary samtools results file"
    self.secondary_sam_file = FileBoxCls(sam_path, "r", fail_msg)
    self.prev_secondary_pair = None
    # get the first pair from the secondary samtools results file
    sam_line = "@"
    try:
      while (sam_line.startswith("@")): #{
        sam_line = self.secondary_sam_file.next()
      #} end while
      self.secondary_pair = P2GPairCls(self.ctg_ids, self.options.read_length,
        self.topology, sam_line, self.log_info)
    except StopIteration, e:
      LogMsg(self,
        "WARNING: no results in secondary samtools results file")
      self.secondary_pair = None
    # end try
    #try:
    #  sam_line = self.secondary_sam_file.next()
    #except StopIteration, e:
    #  LogMsg(self,
    #      "WARNING: no results in secondary samtools results file")
    #  self.secondary_pair = None
    #  return
    ##} end for
    #for sam_line in self.secondary_sam_file: #{
    #  if (not sam_line.startswith("@")): #{
    #    break
    #  #} end if
    ##} end while
    #if (sam_line.startswith("@")): #{
    #  LogMsg(self,
    #      "WARNING: no results in secondary samtools results file")
    #  self.secondary_pair = None
    #else:
    #  self.secondary_pair = P2GPairCls(self.ctg_ids,
    #    self.options.read_length, self.topology, sam_line, self.log_info)
    #} end if
  #} end def

  def ParseSamToolsResults(self, sam_path): #{
    start_time = time.time()
    ExtremeDebugMsg(self, "Parsing samtools results...")
    # set the initial counts
    self.SetInitialCounts()
    # keep track of which pairs have been counted so that
    # no pairs get counted more than once
    self.counted_pairs = dict()
    self.counted_pairs_intron = dict()
    # check each result in the samtools output
    fail_msg = "could not open primary samtools output file"
    for sam_line in FileBoxCls(sam_path, "r", fail_msg): #{
      # skip header lines
      if (sam_line.startswith("@")): #{
        continue
      #} end if
      ExtremeDebugMsg(self, "SAM Line: %s" % sam_line)
      # keep track of the number of reads in the samtools output
      self.num_reads += 1
      pair = P2GPairCls(self.ctg_ids, self.options.read_length, self.topology,
        sam_line, self.log_info)
      # do not count reads more than once and
      # make sure to store the minimum mapq examined for each pair
      if (pair.id in self.counted_pairs): #{
        if (pair.mapq < self.counted_pairs[pair.id]): #{
          self.counted_pairs[pair.id] = pair.mapq
          if (pair.id in self.counted_pairs_intron): #{
            self.counted_pairs_intron[pair.id] = pair.mapq
          #} end if
        #} end if
        ExtremeDebugMsg(self, "Pair already counted: %s MapQ: %i" %
          (pair.id, pair.mapq))
        continue
      #} end if
      ## get the potential supporting pairs from the sam file
      ##potential_pairs = self.GetPotentialPairs(sam_path)
      ## if there are any potentially supporting pairs
      ##if (0 < len(potential_pairs)): #{
      # if the pair is potentially a supporting pair
      if (self.IsPotentialPair(pair)): #{
        ## figure out which pairs are actual supporting pairs
        ##supporting_pairs = self.GetSupportingPairs(potential_pairs)
        ## if there are any actually supporting pairs
        ##if (0 < len(supporting_pairs)): #{
        # if the pair is actually a supporting pair
        if (self.IsSupportingPair(pair)): #{
          self.p2g_all += 1
          ## filter the supporting pairs
          ##self.FilterSupportingPairs(supporting_pairs)
        else:
          ExtremeDebugMsg(self, "--Pair does not support group--")
        #} end if
      else:
        ExtremeDebugMsg(self, "--Pair is not potentially supporting pair--")
      #} end if
    #} end for
    # filter the supporting pairs
    self.FilterSupportingPairs()
    ExtremeDebugMsg(self,
      "Reads: %i, Support: %i (%i intronic), Filtered: %i (%i intronic)\n" %
      (self.num_reads, self.p2g_all, self.p2g_all_intron,
       self.p2g_filt, self.p2g_filt_intron) +
      "Time spent parsing samtools results: %s" % TimeSpent(start_time))
  #} end def

  def SetInitialCounts (self): #{
    self.num_reads       = 0
    self.p2g_all         = 0
    self.p2g_all_intron  = 0
    self.p2g_filt        = 0
    self.p2g_filt_intron = 0
  #} end def

  def IsPotentialPair(self, pair): #{
    ExtremeDebugMsg(self, "Checking whether pair could potentially "
      "support group...")
    # if the mate's chromosome is not right
    if (not pair.CheckMateChromosome(
        self.regions[self.secondary_index].chrom)):
      ExtremeDebugMsg(self, "--wrong chromosome")
      # skip the read
      return False
    #} end if
    # check the read flag to see whether the pair might be supporting
    if (pair.CheckReadFlag()): #{
      ExtremeDebugMsg(self, "Potential read:\n  %s" % pair.ToString())
      return True
    #} end if
    return False
  #} end def

  def GetPotentialPairs(self, sam_path): #{
    ExtremeDebugMsg(self, "Getting potential pairs from samtools results...")
    potential_pairs = dict()
    fail_msg = "could not open samtools output file"
    sam_file = FileBoxCls(sam_path, "r", fail_msg)
    for sam_line in sam_file: #{
      # skip header lines
      if (sam_line.startswith("@")): #{
        continue
      #} end if
      ExtremeDebugMsg(self, "SAM Line: %s" % sam_line)
      self.num_reads += 1
      pair = P2GPairCls(self.ctg_ids, self.options.read_length,
        self.topology, sam_line, self.log_info)
      # if the mate's chromosome is not right
      if (not pair.CheckMateChromosome(
          self.regions[self.secondary_index].chrom)):
        ExtremeDebugMsg(self, "--wrong chromosome")
        # skip the read
        continue
      #} end if
      # check the read flag to see whether the pair might be supporting
      if (pair.CheckReadFlag()): #{
        ExtremeDebugMsg(self, "Potential read:\n  %s" % pair.ToString())
        if (pair.id not in potential_pairs): #{
          potential_pairs[pair.id] = dict()
        #} end if
        potential_pairs[pair.id][pair.num] = pair
      #} end if
    #} end for
    sam_file.close()
    ExtremeDebugMsg(self, "# potential pairs found: %i" % len(potential_pairs))
    return potential_pairs
  #} end def

  def IsSupportingPair(self, pair): #{
    ExtremeDebugMsg(self, "Checking whether pair actually supports group...")
    if ("gap" in self.topology): #{
      return self.CheckGapRegion(pair)
    elif (SAM_BOTH == self.what_to_sam): #{
      # use samtools to get the reads in the secondary region
      return self.CheckRegionWithSam(pair)
    else:
      # just consider the mleft value
      return self.CheckRegionWithMLeft(pair)
    #} end if
  #} end def

  def GetSupportingPairs(self, potential_pairs): #{
    ExtremeDebugMsg(self, "Getting supporting pairs from samtools results...")
    if ("gap" in self.topology): #{
      supporting_pairs = self.CheckGapRegion(potential_pairs)
    elif (SAM_BOTH == self.what_to_sam): #{
      # use samtools to get the reads in the secondary region
      supporting_pairs = self.CheckRegionWithSam(potential_pairs)
    else:
      # looks for reads in the potential pairs hash with
      # the correct mleft value
      supporting_pairs = self.CheckRegionWithMLeft(potential_pairs)
    #} end if
    ExtremeDebugMsg(self, "# supporting pairs found: %i" %
      len(supporting_pairs))
    return supporting_pairs
  #} end def

  def CheckGapRegion(self, pair): #{
    ExtremeDebugMsg(self, "Checking whether pair supports gap group\n"
                   "READ ID: %s" % pair.id)
    region = self.regions[self.secondary_index]
    ExtremeDebugMsg(self,
      "REGION: %i-%i\n" % (region.left, region.right) +
      "%s: mleft = %i" % (pair.id, pair.mleft))
    # get "exonic fragment length"
    try:
      pair.GetExonicFragmentLength(region.blocks)
      ExtremeDebugMsg(self, "exonic fragment length: %i" % pair.frag_len)
    except P2GFragmentError, e:
      raise P2GGroupError("ERROR: could not get exonic fragment length for "
        "ctgs: %s, with blocks: %s\n  %s" % (ctg_ids, region.blocks, e))
    # end try
    if (self.ShouldReadBeCounted(pair, region)): #{
      if (pair.isintronic): #{
        self.p2g_all_intron += 1
        self.counted_pairs_intron[pair.id] = pair.mapq
      #} end if
      self.counted_pairs[pair.id] = pair.mapq
      return True
    #} end if
    return False
  #} end def

  def ShouldReadBeCounted(self, pair, region): #{
    min_len = (1 - self.options.frag_fract) * self.options.frag_len
    max_len = (1 + self.options.frag_fract) * self.options.frag_len
    if ("gap-tandem-duplication" == self.topology): #{
      if (abs(pair.frag_len) <= min_len): #{
        return True
      #} end if
    elif ("gap-nontandem-duplication" == self.topology): #{
      if (abs(pair.frag_len) <= min_len or
          abs(pair.frag_len) >  max_len or
          "RF" == pair.orientation):
        return True
      #} end if
    elif (self.topology in
      ["gap-tandem-inverted_duplication",
       "gap-nontandem-inverted_duplication"]):
      if (abs(pair.frag_len) <= min_len or
          pair.orientation in ["FF", "RR"]):
        return True
      #} end if
    elif ("gap-internal_inversion" == self.topology): #{
      # count read-pairs directed "in" to the inverted region
      if ("RR" == pair.orientation): #{
        if (region.right < pair.MRight()): #{
          return True
        #} end if
      elif ("FF" == pair.orientation): #{
        if (region.left > pair.mleft): #{
          return True
        #} end if
      #} end if
    else:
      LogMsg(self, "WARNING: unrecognized alignment topology: %s" %
        self.topology)
    #} end if
    return False
  #} end def

  def CheckRegionWithSam(self, pair): #{
    ExtremeDebugMsg(self, "Using samtools to check secondary region...")
    # check the previous secondary pair
    if (None != self.prev_secondary_pair): #{
      if (pair.id == self.prev_secondary_pair.id): #{
        ExtremeDebugMsg(self, "Checking previous secondary pair: %s" %
          self.prev_secondary_pair.id)
        try:
          mpair_num = self.prev_secondary_pair.MatePairNum()
          ExtremeDebugMsg(self, "PRIMARY NUM: %s, PREV SECONDARY MATE NUM: "
            "%s" % (pair.num, mpair_num))
        except P2GGroupError, e:
          raise P2GGroupError \
            ("error checking secondary region with samtools for group "
             "%i with ctgs: %s\n  %s" % (self.group_id, self.ctg_ids, e))
        # end try
        # if the current read is the correct mate for the secondary read
        if (pair.num == mpair_num): #{
          # make sure to store the minimum mapq examined for each pair
          self.counted_pairs[pair.id] = min(
            pair.mapq, self.prev_secondary_pair.mapq)
          ExtremeDebugMsg(self, "Counting pair")
          # count the pair
          return True
        #} end if
      else:
        ExtremeDebugMsg(self, "Discarding previous secondary pair: %s" %
          self.prev_secondary_pair.id)
        self.prev_secondary_pair = None
      #} end if
    #} end if
    # look for reads in region B with mates in the potential pairs hash
    if (None == self.secondary_pair): #{
      ExtremeDebugMsg(self, "no (more) reads in secondary sam results file")
      return False
    #} end if
    # skip reads from the secondary file that come
    # before the current primary read
    while (self.secondary_pair.id < pair.id): #{
      ExtremeDebugMsg(self, "Skipping unmatched secondary read: %s" %
        self.secondary_pair.id)
      try:
        self.NextSecondaryPair()
      except StopIteration, e:
        ExtremeDebugMsg(self, "Reached end of secondary sam results file")
        return False
      # end try
    #} end while
    ExtremeDebugMsg(self, "Current secondary pair: %s" %
      self.secondary_pair.id)
    # check all reads from the secondary file that have
    # the same id as the current primary read
    while (self.secondary_pair.id == pair.id): #{
      ExtremeDebugMsg(self, "Checking secondary read: %s (%s)" %
        (self.secondary_pair.id, self.secondary_pair.num))
      try:
        mpair_num = self.secondary_pair.MatePairNum()
        ExtremeDebugMsg(self, "PRIMARY NUM: %s, SECONDARY MATE NUM: %s" %
          (pair.num, mpair_num))
      except P2GGroupError, e:
        raise P2GGroupError \
          ("error checking secondary region with samtools for group "
           "%i with ctgs: %s\n  %s" % (self.group_id, self.ctg_ids, e))
      # end try
      # if the current read is the correct mate for the secondary read
      if (pair.num == mpair_num): #{
        # make sure to store the minimum mapq examined for each pair
        self.counted_pairs[pair.id] = min(pair.mapq, self.secondary_pair.mapq)
        ExtremeDebugMsg(self, "Counting pair")
        self.prev_secondary_pair = None
        # count the pair
        return True
      else:
        self.prev_secondary_pair = self.secondary_pair
        ExtremeDebugMsg(self, "pair numbers do not match, not counting pair")
      #} end if
      try:
        self.NextSecondaryPair()
      except StopIteration, e:
        ExtremeDebugMsg(self, "Reached end of secondary sam results file")
        return False
      # end try
    #} end while
    return False
  #} end def

  def NextSecondaryPair(self): #{
    sam_line = self.secondary_sam_file.next()
    ExtremeDebugMsg(self, "Getting next secondary pair...\n  %s" % sam_line)
    self.secondary_pair = P2GPairCls(self.ctg_ids, self.options.read_length,
      self.topology, sam_line, self.log_info)
    ExtremeDebugMsg(self, "NEW PAIR: %s" % self.secondary_pair.id)
  #} end def

  def CheckRegionWithMLeft(self, pair): #{
    ExtremeDebugMsg(self, "Using mleft to check secondary region...")
    region = self.regions[self.secondary_index]
    ExtremeDebugMsg(self,
      "REGION: %i-%i\n" % (region.left, region.right) +
      "%s: mleft = %i" % (pair.id, pair.mleft))
    if (region.left <= pair.mleft and pair.mleft <= region.right): #{
      ExtremeDebugMsg(self, "Counting pair")
      # make sure to store the minimum mapq examined for each pair
      self.counted_pairs[pair.id] = pair.mapq
      # count the pair
      return True
    #} end if
    return False
  #} end def

  def FilterSupportingPairs(self): #{
    for read_id, mapq in self.counted_pairs.iteritems(): #{
      ExtremeDebugMsg(self, "%s mapq: %i" % (read_id, mapq))
      if (self.options.min_mapq < mapq): #{
        self.p2g_filt += 1
      #} end if
    #} end for
    for read_id, mapq in self.counted_pairs_intron.iteritems(): #{
      ExtremeDebugMsg(self, "%s mapq: %i intronic" % (read_id, mapq))
      if (self.options.min_mapq < mapq): #{
        self.p2g_filt_intron += 1
      #} end if
    #} end for
  #} end def

  def ToString(self): #{
    data_fields = list([
      "%i" % self.group_id,
      self.coords[UPSTREAM_INDEX].ToString(),
      self.coords[DOWNSTREAM_INDEX].ToString(),
      self.topology,
      str(self.ctg_overlap),
      self.ctg_ids,
    ])
    return "\t".join(data_fields)
  #} end def

  def ParseSupportString(self, support_string): #{
    (self.group_id, self.num_reads,
     self.p2g_exonic, self.p2g_all,
     self.p2g_filt_exonic, self.p2g_filt_all) = \
      map (int, support_string.split("\t"))
  #} end if

  def SupportString(self): #{
    data_fields = list([
      self.group_id,
      self.num_reads,
      self.p2g_all - self.p2g_all_intron,
      self.p2g_all,
      self.p2g_filt - self.p2g_filt_intron,
      self.p2g_filt,
    ])
    return "\t".join(map(str, data_fields))
  #} end def
#} end class

class PairToGenomeCoordsCls: #{
  def __init__(self, coords_string): #{
    self.ParseCoordsString(coords_string)
  #} end def

  def ParseCoordsString(self, coords_string): #{
    (chrom, left_str, right_str, self.breakside,
      self.blocks) = coords_string.split(" ")
    self.chrom = NormalizeChrID(chrom)
    self.left  = int(left_str)
    self.right = int(right_str)
  #} end def

  def Span(self): #{
    return abs(self.left - self.right) + 1
  #} end def

  def ToString(self, use_chr=False): #{
    if (use_chr): #{
      chromosome = AddChr(self.chrom)
    else:
      chromosome = self.chrom
    #} end if
    data_fields = list([
      chromosome,
      "%i" % self.left,
      "%i" % self.right,
      self.breakside,
      self.blocks,
    ])
    return " ".join(data_fields)
  #} end def
#} end class

class SearchRegionCls: #{
  def __init__(self, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    self.left  = -1
    self.right = -1
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def GetCoords(self, topology, full_coords, length): #{
    # set chromosome & breakpoint side
    self.chrom = full_coords.chrom
    self.breakside = full_coords.breakside
    # use the whole region for gap groups
    if ("gap" in topology): #{
      self.left  = full_coords.left
      self.right = full_coords.right
    else:
      self.length = length
      # if the whole region is shorter than the desired length
      if (full_coords.Span() < length): #{
        # extend the full coords to the correct length
        self.ExtendFullRegion(full_coords)
      else:
        # walk along the blocks to get the region
        self.GetCoordsByWalk(full_coords)
      #} end if
    #} end if
    # ensure that left side is not negative
    self.left = max(0, self.left)
  #} end def

  def ExtendFullRegion(self, full_coords): #{
    ExtremeDebugMsg(self, "Extending full coordinates...");
    if ("left" == self.breakside): #{
      self.left  = full_coords.left
      self.right = full_coords.left + self.length
    else: # "right" == self.breakside
      self.left  = full_coords.right - self.length
      self.right = full_coords.right
    #} end if
  #} end def

  def GetCoordsByWalk(self, full_coords): #{
    ExtremeDebugMsg(self,
      "Walking from %s...\n" % self.breakside +
      "BLOCKS: %s" % full_coords.blocks)
    if ("left" == self.breakside): #{
      self.left  = full_coords.left
      self.right = self.WalkForwards(full_coords.left, full_coords.blocks)
    else: # "right" == self.breakside
      self.left  = self.WalkBackwards(full_coords.right, full_coords.blocks)
      self.right = full_coords.right
    #} end if
  #} end def

  def WalkForwards(self, walk_start, blocks_string): #{
    walk_len = self.length
    walk_end = -1
    blocks = blocks_string.split(",")
    if (1 == len(blocks)): #{
      (block_left, block_right) = map(int, blocks[0].split("-"))
      walk_end = block_left + walk_len
    else:
      for block in blocks: #{
        (block_left, block_right) = map(int, block.split("-"))
        block_len = abs(block_left-block_right) + 1
        # walk over the current block
        if (block_len < walk_len): #{
          walk_len -= block_len
        # end on the current block
        else:
          walk_end = block_left + walk_len
          walk_len = 0
          break
        #} end if
      #} end for
      # if there is still some walking to do
      if (0 < walk_len): #{
        # use the right side of the last block
        walk_end = int(blocks[-1].split("-")[1])
      #} end if
    #} end if
    if (0 > walk_end): #{
      raise P2GGroupError("Could not get walk end with blocks: %s" %
        blocks_string)
    #} end if
    return walk_end
  #} end def

  def WalkBackwards(self, walk_start, blocks_string): #{
    walk_len = self.length
    walk_end = -1
    blocks = blocks_string.split(",")
    if (1 == len(blocks)): #{
      (block_left, block_right) = map(int, blocks[0].split("-"))
      walk_end = block_right - walk_len
    else:
      for block in reversed(blocks): #{
        (block_left, block_right) = map(int, block.split("-"))
        block_len = abs(block_left-block_right) + 1
        if (block_len < walk_len): #{
          walk_len -= block_len
        else:
          walk_end = block_right - walk_len
          walk_len = 0
          break
        #} end if
      #} end for
      # if there is still some walking to do
      if (0 < walk_len): #{
        # use the left side of the first block
        walk_end = int(blocks[0].split("-")[0])
      #} end if
    #} end if
    if (0 > walk_end): #{
      raise P2GGroupError("Could not get walk end with blocks: %s" %
        blocks_string)
    #} end if
    return walk_end
  #} end def

  def ToString(self, use_chr=False): #{
    if (use_chr): #{
      chromosome = AddChr(self.chrom)
    else:
      chromosome = self.chrom
    #} end if
    return "%s:%i-%i" % (chromosome, self.left, self.right)
  #} end def
#} end class

class P2GPairCls: #{
  def __init__(self, ctg_ids, read_length, topology, sam_line,
      log_info=None):
    self.log_info    = log_info
    self.read_length = read_length
    self.topology    = topology
    self.isintronic  = False
    self.ParseSAMLine(sam_line)
    self.SetPairNum()
  #} end def

  def ParseSAMLine(self, sam_line): #{
    (self.id, flag, chromosome, left, mapq, self.cigar, self.mchrom, mleft,
      genomic_frag_len) = sam_line.split("\t")[:9]
    self.chrom = NormalizeChrID(chromosome)
    try:
      self.flag  = int(flag)
      self.left  = int(left)
      self.mapq  = int(mapq)
      self.mleft = int(mleft)
      self.genomic_frag_len = int(genomic_frag_len)
    except ValueError, e:
      raise P2GGroupError("error parsing SAM line:\n%s\n%s" % (sam_line, e))
    #} end if
  #} end def

  def SetPairNum(self): #{
    if (self.flag & 0x0040): #{
      self.num = "first"
    # if the read is the second read in the pair
    elif (self.flag & 0x0080): #{
      self.num = "second"
    # otherwise
    else:
      LogMsg(self,
        "WARNING: read was neither the first nor second pair: %s" % self.id)
      self.num = None
    #} end if
  #} end def

  def Right(self): #{
    return (self.left + self.read_length)
  #} end def

  def MRight(self): #{
    return (self.mleft + self.read_length)
  #} end def

  def DownStreamLeft(self): #{
    return min(self.left, self.mleft)
  #} end def

  def UpStreamLeft(self): #{
    return max(self.left, self.mleft)
  #} end def

  def MatePairNum(self): #{
    if (not hasattr(self, "num") or
        None == self.num):
      # set the read num first
      self.SetPairNum()
    #} end if
    # if the read is the first read in the pair
    if ("first" == self.num): #{
      # its mate is the second
      return "second"
    # if the read is the second read in the pair
    elif ("second" == self.num): #{
      # its mate is the first
      return "first"
    # otherwise
    else:
      raise P2GGroupError \
        ("read was neither the first nor second pair: %s, " % self.id +
         "there may be a problem with the read alignment file.")
    #} end if
  #} end def

  def CheckMateChromosome(self, desired_chrom): #{
    ExtremeDebugMsg(self, "Mine: %s, Desired: %s, Mate: %s" %
      (self.chrom, desired_chrom, self.mchrom))
    if (desired_chrom == self.chrom): #{
      desired_chrom = "="
    #} end if
    return (self.mchrom == desired_chrom)
  #} end def

  def CheckReadFlag(self): #{
    ExtremeDebugMsg(self, "Checking flag of read...")
    # check that the paired-read flag is properly set
    if (self.flag & 0x0001): #{
      # if the read or its mate is unmapped
      if (self.flag & 0x0004 or self.flag & 0x0008): #{
        # skip the read
        ExtremeDebugMsg(self, "--read or pair is unmapped: %s" % self.id)
        return False
      #} end if
      # strand 0 is the forward strand, strand 1 is reverse
      # assume both reads are on the forward strand
      (strand, mstrand) = (0,0)
      if (self.flag & 0x0010): #{
        strand = 1
      #} end if
      if (self.flag & 0x0020): #{
        mstrand = 1
      #} end if
      ExtremeDebugMsg(self, "Strand: %i, MStrand: %i" % (strand, mstrand))
      reads = list()
      # if the read is the first read in the pair
      if ("first" == self.num): #{
        self.num = "first"
        reads.append(P2GReadCls(strand,  self.left))
        reads.append(P2GReadCls(mstrand, self.mleft))
      # if the read is the second read in the pair
      elif ("second" == self.num): #{
        self.num = "second"
        reads.append(P2GReadCls(mstrand, self.mleft))
        reads.append(P2GReadCls(strand,  self.left))
      # otherwise
      else:
        LogMsg(self,
          "WARNING: read was neither the first nor second pair: %s" % self.id)
        return False
      #} end if
      if (self.CheckStrandAndOrientation(reads)): #{
        #LogMsg(self, "Good Pair: %s" % self.id)
        return True
      #} end if
    else:
      LogMsg(self,
        "WARNING: read was not paired in sequencing: %s" % self.id)
    #} end if
    return False
  #} end def

  def CheckStrandAndOrientation(self, reads): #{
    ExtremeDebugMsg(self, "Checking strand and orientation of read...")
    if ("interchr" == self.topology): #{
      self.orientation = "NA"
      # count every pair
      return True
    #} end if
    # determine the read-pair orientation
    self.GetOrientation(reads)
    if (self.topology in ["local-inversion", "intrachr-opp-strand"]): #{
      # count read pairs with the same strand
      if (reads[0].strand == reads[1].strand): #{
        return True
      #} end if
    elif (self.topology in
      ["junction-duplication", "intrachr-non-colinear", "end-duplication"]):
      # count read pairs with outwards (RF) orientation
      if ("RF" == self.orientation): #{
        return True
      #} end if
    elif (self.topology in
      ["read-through", "intrachr-same-strand", "gap-nontandem-duplication"]):
      # count read pairs with opposite strands
      if (reads[0].strand != reads[1].strand): #{
        return True
      #} end if
    elif ("gap-tandem-duplication" == self.topology): #{
      # count read pairs with inwards (FR) orientation
      if ("FR" == self.orientation): #{
        return True
      #} end if
    elif (self.topology in
      ["gap-tandem-inverted_duplication",
       "gap-nontandem-inverted_duplication"]):
      # count pairs that do not have outwards (RF) orientation
      if ("RF" != self.orientation): #{
        return True
      #} end if
    elif ("gap-internal_inversion" == self.topology): #{
      # count left- or right-spooned reads
      if (self.orientation in ["RR", "FF"]): #{
        return True
      #} end if
    elif ("gap-genic_rearrangement" == self.topology): #{
      # REMINDER: how to deal with gap genic-rearrangments?
      return True # JUST FOR NOW
    else:
      LogMsg(self, "WARNING: in CheckStrandAndOrientation: "
        "contigs %s have unrecognized topology: %s" %
        (self.ctg_ids, self.topology))
    #} end if
    return False
  #} end def

  def GetOrientation(self, reads): #{
    if (reads[0].strand == reads[1].strand): #{
      # if the first read is on the forward strand
      if (0 == reads[0].strand): #{
        self.orientation = "FF"
      # if the first read is on the reverse strand
      else: # 1 == reads[0].strand
        self.orientation = "RR"
      #} end if
    else: # reads[0].strand != reads[1].strand
      # assume the reads are oriented inwards
      self.orientation = "FR"
      # if the first read is on the forward strand
      if (0 == reads[0].strand): #{
        # and the second read is upstream
        if (reads[0].left - reads[1].left + self.read_length > 0): #{
          # the reads are oriented outwards
          self.orientation = "RF"
        #} end if
      # if the first read is on the reverse strand
      else: # 1 == reads[0].strand
        # and the second read is downstream
        if (reads[1].left - reads[0].left + self.read_length > 0): #{
          # the reads are oriented outwards
          self.orientation = "RF"
        #} end if
      #} end if
    #} end if
    ExtremeDebugMsg(self, "Read orientation: %s" % self.orientation)
  #} end def

  def GetExonicFragmentLength(self, blocks_str): #{
    fragment = P2GFragmentCls(self, log_info=self.log_info)
    ExtremeDebugMsg(self,
      "Fragment: %s\n" % fragment.LeftsString() +
      "BLOCKS: %s" % blocks_str)
    for block in blocks_str.split(","): #{
      (block_left, block_right) = map(int, block.split("-"))
      ExtremeDebugMsg(self,
        "BLOCK: %s; Left: %i, Right: %i" % (block, block_left, block_right))
      # skip blocks occurring before the fragment
      if (fragment.IsBlockUpstream(block_right)): #{
        continue
      #} end if
      # skip blocks occurring after the fragment
      if (fragment.IsBlockDownstream(block_left)): #{
        fragment.EndBetweenBlocks()
        break
      #} end if
      # make sure to start counting at the left side of the fragment
      if (fragment.IsBlockFirst(block_left)): #{
        block_left = fragment.downstream_left
      #} end if
      # if the fragment ends before the block,
      # add only to the end of the fragment
      if (fragment.IsBlockLast(block_right)): #{
        fragment.EndInBlock(block_left)
        break
      # if the block ends before the fragment,
      # add the whole block length
      else: # block_right < self.upstream_left
        fragment.AddBlock(block_left, block_right)
      #} end if
      fragment.prev_right = block_right
    #} end for
    fragment.CheckIfBlocksMiss()
    fragment.CheckLength()
    self.isintronic = fragment.isintronic
    self.frag_len = fragment.Length()
  #} end def

  def ToString(self): #{
    return " ".join([
      "%s(%s)"  % (self.id, self.num),
      "%s:%i"   % (self.chrom, self.left),
      "%s:%i"   % (self.mchrom, self.mleft),
      "%s"      % self.orientation,
      "mapq:%i" % self.mapq,
    ])
  #} end if
#} end class

class P2GReadCls: #{
  def __init__(self, strand, left_coord): #{
    self.strand = strand
    self.left   = left_coord
  #} end def
#} end class

class P2GFragmentCls: #{
  def __init__(self, pair, log_info=None): #{
    self.log_info    = log_info
    self.read_length = pair.read_length
    self.downstream_left = pair.DownStreamLeft()
    self.upstream_left   = pair.UpStreamLeft()
    self.frag_len = 0
    # assume that the first read in the pair is aligned between blocks
    self.isintronic = True
    self.prev_right = self.downstream_left
  #} end def

  def LeftsString(self): #{
    return "%i-%i" % (self.downstream_left, self.upstream_left)
  #} end def

  def IsBlockUpstream(self, block_right): #{
    return (block_right < self.downstream_left)
  #} end def

  def IsBlockDownstream(self, block_left): #{
    return (self.upstream_left < block_left)
  #} end def

  def IsBlockFirst(self, block_left): #{
    if (0 == self.frag_len): #{
      # check whether the first read in the pair is aligned between blocks
      if (block_left < self.downstream_left): #{
        self.isintronic = False
      #} end if
      return True
    #} end if
    return False
  #} end def

  def IsBlockLast(self, block_right): #{
    return (block_right >= self.upstream_left)
  #} end def

  def EndBetweenBlocks(self): #{
    self.AddBlock(self.prev_right, self.upstream_left)
    # the second read in the pair is aligned between blocks
    self.isintronic = True
  #} end def

  def EndInBlock(self, block_left): #{
    self.AddBlock(block_left, self.upstream_left)
  #} end def

  def AddBlock(self, block_left, block_right): #{
    block_len = (block_right - block_left) + 1
    self.frag_len += block_len
    ExtremeDebugMsg(self, "Adding %i, length = %i" % (block_len, self.frag_len))
  #} end def

  def CheckIfBlocksMiss(self): #{
    if (self.prev_right <= self.downstream_left and 0 == self.frag_len): #{
      self.frag_len = (self.upstream_left - self.downstream_left) + 1
    #} end if
  #} end def

  def CheckLength(self): #{
    if (1 > self.frag_len): #{
      raise P2GFragmentError("Could not get exonic fragment length for region: "
          "%i-%i " % (self.downstream_left, self.upstream_left))
    #} end if
  #} end def

  def Length(self): #{
    return (self.frag_len + self.read_length)
  #} end def
#} end class


#### EXCEPTION CLASSES ####
class P2GGroupError(MyError): #{
  pass
#} end class

class P2GFragmentError(MyError): #{
  pass
#} end class
