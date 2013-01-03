#! /usr/bin/env python
"""
realigner.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os, time

# import custom modules
from utils.error import MyError
from utils.general import TimeSpent, GetCommand
from utils.messages import LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import FileBoxCls
from utils.subprocesses import RunCommandFromList
from alignment_processing.alignment_functions import (ParseAlignmentFile,
  CalcOverlap, FixNonChromAlign, NoAlignmentsError)
from parsers.fasta import FastaFileCls
from common.coord_pair import CoordPairCls

# constants
MIN_REGION_OVERLAP_FRACTION = 0.70
MAX_MATCH_DIFF = 10

class RealignerCls: #{
  def __init__(self, options, log_info): #{
    self.options = options
    self.log_info = log_info
    # contigs[ctg_id] = realign result
    #   realign result: best alignment of contig to transcripts
    self.contigs = dict()
    self.missing = set()
    #self.contig_seqs = dict()
    self.SetupPaths()
  #} end def

  def SetupPaths(self): #{
    query_file_name = "%s.event_contigs.fa" % self.options.lib
    self.query_path = os.path.join(self.options.realign_dir,
      query_file_name)
    results_file_name = "%s.realign.psl" % self.options.lib
    self.results_path = os.path.join(self.options.realign_dir,
      results_file_name)
    fine_results_file_name = "%s.realign.fine.psl" % self.options.lib
    self.fine_results_path = os.path.join(self.options.realign_dir,
      fine_results_file_name)
  #} end def

  #def UpdateContigs(self, group, good_members, store_seq): #{}
  def UpdateContigs(self, group, good_members, event_type): #{
    for index in good_members: #{
      ctg_id = group.members[index].contig_info.id
      if (ctg_id not in self.contigs): #{
        self.contigs[ctg_id] = RealignedContigCls(ctg_id,
          log_info=self.log_info)
      #} end if
      self.contigs[ctg_id].UpdateMembers(group.members[index], event_type)
      self.missing.add(ctg_id)
      #if (store_seq): #{
      #  self.contig_seqs[ctg_id] = None
      #} end if
    #} end for
  #} end def

  def RealignContigs(self): #{
    start = time.time()
    if (not os.path.isfile(self.results_path)): #{
      LogMsg(self, "Realignment file does not exist: \"%s\"" %
        self.results_path)
      self.options.use_existing_realigns = False
    #} end if
    if (not os.path.isfile(self.fine_results_path)): #{
      LogMsg(self, "Fine realignment file does not exist: \"%s\"" %
        self.fine_results_path)
      self.options.use_existing_realigns = False
    #} end if
    if (self.options.use_existing_realigns): #{
      LogMsg(self, "Using existing contig-to-transcript realignment files.")
      self.LoadContigSequences()
    else:
      LogMsg(self, "Performing contig realignment...")
      # create a query file containing the correct contig sequences
      self.CreateQueryFile()
      # run the blat alignment without the "-fine" option
      self.RunBlatAlignment(self.results_path, use_fine=False)
      # run the blat alignment with the "-fine" option
      self.RunBlatAlignment(self.fine_results_path, use_fine=True)
    #} end if
    # parse the alignment results
    self.ParseRealignments(self.results_path)
    self.ParseRealignments(self.fine_results_path)
    LogMsg(self, "Time spent processing contig realignments: %s" %
      TimeSpent(start))
  #} end def

  def LoadContigSequences(self): #{
    LogMsg(self, "Loading contig sequences from %s..." % self.query_path)
    contigs_file = FastaFileCls(self.query_path)
    num_loaded = 0
    for contig in contigs_file: #{
      if (contig.id in self.contigs): #{
        self.contigs[contig.id].sequence = contig.sequence.lower()
        num_loaded += 1
      #} end if
    #} end for
    if (num_loaded != len(self.contigs)): #{
      LogMsg(self, "WARNING: only loaded %i of %i contig sequences! " %
        (num_loaded, len(self.contigs)))
      #  ",".join(missed))
    #} end if
  #} end def

  def CreateQueryFile(self): #{
    LogMsg(self, "Creating query file...")
    query_file = FileBoxCls(self.query_path, "w", "cannot create query "
      "contig sequences file")
    all_contigs_file = FileBoxCls(self.options.ctg_seq_path, "r",
      "cannot read contig sequences file")
    seqs_found = False
    num_written = 0
    for id_line in all_contigs_file: #{
      seq_line = all_contigs_file.next()
      if (not id_line.startswith(">")): #{
        raise RealignerError("invalid contig id line in sequece file:\n%s" %
          id_line)
      #} end if
      # extract the contig id from the line
      ctg_id = id_line.lstrip(">").split()[0]
      #DebugMsg(self, "Contig ID from sequence file: %s" % ctg_id)
      # if the contig is represented in one of the potential predictions
      if (ctg_id in self.contigs): #{
        #DebugMsg(self, "Writing sequence to query file.")
        # write it to the query file
        query_file.WriteLine(id_line)
        query_file.WriteLine(seq_line)
        seqs_found = True
        num_written += 1
        self.contigs[ctg_id].written = True
        self.contigs[ctg_id].sequence = seq_line.lower()
        #if ("itd" in self.contigs[ctg_id].types): #{
        #  self.contigs[ctg_id].sequence = seq_line.lower()
        #} end if
        self.missing.discard(ctg_id)
      #} end if
      #if (ctg_id in self.contig_seqs): #{
      #  self.contig_seqs[ctg_id] = seq_line
      #} end if
    #} for
    if (not seqs_found): #{
      raise RealignerError("could not find any contig sequences in %s" %
        self.options.ctg_seq_path)
    #} end if
    if (num_written != len(self.contigs)): #{
      #missed = list()
      #for contig in self.contigs.itervalues(): #{
      #  if (not contig.written): #{
      #    missed.append(contig.id)
      #  #} end if
      #} end for
      LogMsg(self, "WARNING: only wrote %i of %i contig sequences! " %
        (num_written, len(self.contigs)) + "Missing: %s" %
        ",".join(sorted(self.missing)))
      #  ",".join(missed))
    #} end if
    all_contigs_file.Close()
    query_file.Close()
  #} end def

  def RunBlatAlignment(self, results_path, use_fine=False): #{
    LogMsg(self, "Realigning contigs to transcript sequences...")
    start = time.time()
    # use "stepSize=5" and "minScore=4" to allow short matches
    # use "repMatch=2253" to mimic webBLAT results
    command = [GetCommand(self, "blat"), "-stepSize=5", "-repMatch=2253",
      "-minScore=4", self.options.tran_seq_path, self.query_path, results_path]
    # use "fine" option to avoid strange misalignments around single mismatches
    if (use_fine): #{
      command.append("-fine")
    #} end if
    DebugMsg(self, " ".join(command))
    return_code = RunCommandFromList(command, dpt=self.options.dpt)
    if (0 != return_code and not self.log_info['debug']): #{
      LogMsg(self, " ".join(command))
    #} end if
    if (0 > return_code): #{
      raise RealignerError("Realignment command was terminated by signal "
        "%i" % return_code)
    elif (0 < return_code):
      raise RealignerError("Error running realignment command: %i" %
        return_code)
    #} end if
    LogMsg(self, "Time spent performing realignment: %s" %
      TimeSpent(start))
  #} end def

  def ParseRealignments(self, results_path): #{
    LogMsg(self, "Parsing realignment results...")
    #filters = { 'count': 75, 'identity': 70.0 }
    filters = { 'identity': 40.0 }
    try: #{
      aligns = ParseAlignmentFile(results_path, filters, self.log_info)
    except NoAlignmentsError, e:
      LogMsg(self, "WARNING: no contig-to-transcript alignments were found.")
      return
    #} end try
    for align in map(FixNonChromAlign, aligns): #{
      ctg_id = align.query
      #if (None != self.contigs[ctg_id]): #{
      #  raise RealignerError("multiple alignments returned for contig: %s" %
      #    ctg_id)
      #} end if
      #new_align = MinimalAlignmentCls(align)
      #if (None == self.contigs[ctg_id] or
      #    self.contigs[ctg_id].match < new_align.match): #{
      #  if (None != self.contigs[ctg_id]): #{
      #    DebugMsg(self, "Replacing alignment")
      #  #} end if
      #  self.contigs[ctg_id] = new_align
      #} end if
      self.UpdateContigAligns(self.contigs[ctg_id], align)
    #} end for
  #} end def

  def UpdateContigAligns(self, contig, full_align): #{
    new_align = MinimalAlignmentCls(full_align)
    #DebugMsg(self, "ALIGN %s" % new_align.ToString())
    contig.ctg_length = new_align.ctg_length
    if ("itd" in contig.types): #{
      if (5 > new_align.ctg_start): #{
        contig.start_aligns.append(new_align)
      #} end if
      if (5 > (new_align.ctg_length - new_align.ctg_end)): #{
        contig.end_aligns.append(new_align)
      #} end if
      # perhaps also check that 0 = num_query_gaps?
      #if (self.options.read_length <= new_align.Span()): #{}
      if (self.options.read_length <= new_align.match): #{
        self.CheckForFullDupAligns(new_align, contig)
      #} end if
    #} end if
    if (None == contig.best_align or
        contig.best_align.match < new_align.match): #{
      if (None != contig.best_align): #{
        ExtremeDebugMsg(self, "Replacing alignment (old: %i, new: %i)" %
          (contig.best_align.match, new_align.match))
      #} end if
      contig.best_align = new_align
    #} end if
    if (0 < len(contig.region_pairs)): #{
      for region_pair in contig.region_pairs.itervalues(): #{
        region_pair.UpdateGenes(new_align)
      #} end for
    #} end if
  #} end def

  def CheckForFullDupAligns(self, new_align, contig): #{
    #ExtremeDebugMsg(self, "Checking for full dup align: %s" % new_align.ctg_coords)
    for (member_id, region) in contig.full_dup_regions.iteritems(): #{
      #ExtremeDebugMsg(self, "  Member full dup coords: %s" % region)
      if (new_align.ctg_coords.Contains(region)): #{
        #ExtremeDebugMsg(self, "    Adding full dup target: %s" %
        #    new_align.target)
        if (member_id not in contig.full_dup_aligns): #{
          contig.full_dup_aligns[member_id] = set()
        #} end if
        contig.full_dup_aligns[member_id].add(new_align.target)
      #} end if
    #} end for
    #if (new_align.ExactMatch()): #{
    #  for (member_id, region) in contig.full_dup_regions.iteritems(): #{
    #    if (new_align.ctg_coords.Contains(region)): #{
    #      contig.full_dup_aligns[member_id] = new_align.target
    #    #} end if
    #  #} end for
    #} end if
  #} end def
#} end class

class RealignedContigCls: #{
  def __init__(self, ctg_id, log_info=None): #{
    self.id = ctg_id
    self.types = set()
    # region_pairs[member_id] = RegionPairCls
    self.region_pairs = dict()
    # full_dup_regions[member_id] = CoordPairCls
    #  (union of gap-coords and dup-coords)
    self.full_dup_regions = dict()
    # full_dup_aligns[member_id] = set of genes to which contig has
    #  alignment covering the "full_dup_region" for that member
    self.full_dup_aligns = dict()
    self.sequence = None
    self.best_align = None
    self.start_aligns = list()
    self.end_aligns = list()
    self.ctg_length = None
    self.written = False
    self.log_info = log_info
  #} end def

  def UpdateMembers(self, member, event_type): #{
    self.types.add(event_type)
    if ("fusion" == event_type): #{
      # update the region pairs
      self.region_pairs[member.IDString()] = RegionPairCls(member,
        log_info=self.log_info)
    elif ("itd" == event_type and member.gap and member.GapIsInternal()):
      # add a full_dup_region
      full_dup_region = CoordPairCls(member.align_info_B.ctg_coords)
      full_dup_region.Union(member.meta_fields['gap_coords'])
      full_dup_region.Union(member.meta_fields['dup_coords'])
      self.full_dup_regions[member.IDString()] = full_dup_region
    #} end if
  #} end def
#} end class

class RegionPairCls: #{
  def __init__(self, member, log_info=None): #{
    self.member_id = member.IDString()
    # extend regions to beginning and end of contig, respectively
    #self.starts = {'A': member.align_info_A.ctg_start,
    self.starts = {'A': 1, 'B': member.align_info_B.ctg_start}
    self.ends   = {'A': member.align_info_A.ctg_end,
      'B': member.contig_info.length}
      #'B': member.align_info_B.ctg_end}
    self.best_aligns = {'A': None, 'B': None}
    self.best_matches = {'A': 0, 'B': 0}
    self.gene_sets = {'A': set(), 'B': set()}
    self.log_info = log_info
  #} end def

  def UpdateGenes(self, new_align): #{
    ExtremeDebugMsg(self, "Updating genes for member %s regions" %
      self.member_id)
    #align_span = (new_align.ctg_end - new_align.ctg_start) + 1
    for region_id in ["A", "B"]: #{
      if (new_align.target in self.gene_sets[region_id] and
          new_align.match < self.best_matches[region_id]): #{
        ExtremeDebugMsg(self, "  Region %s: do not need to process align: "
          "%i-%i, Target: %s." % (region_id, new_align.ctg_start,
          new_align.ctg_end, new_align.target))
        continue
      #} end if
      #min_overlap = (min(self.Span(region_id), align_span) *
      #  MIN_REGION_OVERLAP_FRACTION)
      (overlap, fraction) = CalcOverlap(self.starts[region_id],
        self.ends[region_id], new_align.ctg_start, new_align.ctg_end)
      match_diff = self.best_matches[region_id] - new_align.match
      ExtremeDebugMsg(self, "  Region %s: %i-%i. New align: %i-%i, "
        "Target: %s.\n  Overlap: %i (%.2f), Match Difference: %i" % (region_id,
        self.starts[region_id], self.ends[region_id], new_align.ctg_start,
        new_align.ctg_end, new_align.target, overlap, fraction, match_diff))
      if (MIN_REGION_OVERLAP_FRACTION < fraction and
          MAX_MATCH_DIFF > match_diff): #{
        ExtremeDebugMsg(self, "  Adding target to region %s" % region_id)
        self.gene_sets[region_id].add(new_align.target)
        if (self.best_matches[region_id] < new_align.match): #{
          self.best_matches[region_id] = new_align.match
          self.best_aligns[region_id] = new_align
        #} end if
      #} end if
    #} end for
  #} end def

  def GeneSetsOverlap(self): #{
    if (self.gene_sets['A'].isdisjoint(self.gene_sets['B'])): #{
      return False
    #} end if
    return True
  #} end def
#} end class

# to minimize the space required to store all the alignments, only
# store the necessary information
class MinimalAlignmentCls: #{
  def __init__(self, align): #{
    self.match = align.match
    self.ctg_start  = align.qstart
    self.ctg_end    = align.qend
    self.ctg_coords = CoordPairCls(self.ctg_start, self.ctg_end)
    self.ctg_length = align.query_len
    self.strand     = align.query_strand
    self.target     = align.target
    self.num_query_gaps = align.qnuminsert
    # only save the blocks if there are query gaps that might need checking
    self.query_blocks = None
    if (0 < self.num_query_gaps): #{
      self.query_blocks = align.query_blocks
      # if aligned to the negative strand, need to reverse all the blocks
      if ("-" == align.query_strand): #{
        self.query_blocks.reverse()
        for index in range(len(self.query_blocks)): #{
          self.query_blocks[index].reverse()
        #} end for
      #} end if
    #} end if
    # if no more than 4 bases are not matches, mark a perfect alignment
    if (5 > (self.ctg_length - self.match)): #{
      self.perfect = True
    else:
      self.perfect = False
    #} end if
  #} end def

  def ExactMatch(self): #{
    if (self.match == self.Span()): #{
      return True
    #} end if
    return False
  #} end def

  def Span(self): #{
    return abs(self.ctg_end - self.ctg_start)+1
  #} end def

  def ToString(self): #{
    data_str = "%i-%i(%i) %s(%s)" % (self.ctg_start, self.ctg_end, self.match,
      self.target, self.strand)
    if (0 < self.num_query_gaps): #{
      data_str += ": %s" % ",".join(["%i-%i" % (block[0],block[1]) for
        block in self.query_blocks])
    #} end if
    return data_str
  #} end if
#} end class

#### EXCEPTION CLASSES ####
class RealignerError(MyError): #{
  pass
#} end class
