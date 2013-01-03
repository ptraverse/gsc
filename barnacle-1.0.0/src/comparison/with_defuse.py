#! /usr/bin/env python
"""
with_defuse.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, re, sys, time, traceback

# import custom modules
from version import VERSION
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  CheckConfigCommands, GetCommand)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.multi_dict import (SetupMultiDict, CheckMultiDict, AccessMultiDict,
  TraverseMultiDict, MultiDictError)
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls)
from utils.subprocesses import RunCommandFromList
from common.breakpoint import SortBreakpoints, BreakpointPairsAreEqual
from parsers.candidate_group_parser import CandidateGroupParserCls
from alignment_processing.alignment_functions import (ParseAlignmentFile,
  FixAlign, AlignBlockCls, NoAlignmentsError)
from defuse_prediction import (DEFUSE_COMPARISON_HEADER, ALIGN_RESULTS_HEADER,
  DefusePredictionCls, DefusePredictionError)
from breakpoint_pair import BreakpointPairCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"
COUNTS_HEADER1 = "\t".join(["        ","total      ","lcommon    ",
  "gcommon    ","bcommon    ","total   ","lcommon ","gcommon ","bcommon "])
COUNTS_HEADER2 = "\t".join(["Category","predictions","predictions",
  "predictions","predictions","partners","partners","partners","partners"])
CONTINUOUS_THRESHOLD = 5

class DefuseComparisonCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    if (self.options.align): #{
      CheckConfigCommands(self, "blat")
    #} end if
    # lcom: common by general location (big_buffer)
    # gcom: common by gene partner names
    # bcom: common by specific breakpoint (small_buffer)
    # qual0: total
    # qual1: not mitochondrial
    # qual2: not mitochondrial, not a readthrough-like event
    # qual3: not mitochondrial, not readthrough, and medium or better
    #         confidence (based on min value for 'probability' field)
    # qual4: not mitochondrial, not readthrough, and high confidence
    #         (based on min value for 'probability' field)
    # qual5: not mitochondrial, not readthrough, high confidence, and
    #         good coverage
    self.counts = {
      "defuse": {
        "total": {"qual0":0,    "qual1":0,    "qual2":0,
                  "qual3":0,    "qual4":0,    "qual5":0,},
        "lcom":  {"qual0":0,    "qual1":0,    "qual2":0,
                  "qual3":0,    "qual4":0,    "qual5":0,},
        "gcom":  {"qual0":0,    "qual1":0,    "qual2":0,
                  "qual3":0,    "qual4":0,    "qual5":0,},
        "bcom":  {"qual0":0,    "qual1":0,    "qual2":0,
                  "qual3":0,    "qual4":0,    "qual5":0,},
      },
      "barnacle": {
        "total":0, "lcom":0, "gcom":0, "bcom":0,
      }
    }
    self.partners = {
      "defuse": {
        "total": {"qual0":set(),"qual1":set(),"qual2":set(),
                  "qual3":set(),"qual4":set(),"qual5":set(),},
        "lcom":  {"qual0":set(),"qual1":set(),"qual2":set(),
                  "qual3":set(),"qual4":set(),"qual5":set(),},
        "gcom":  {"qual0":set(),"qual1":set(),"qual2":set(),
                  "qual3":set(),"qual4":set(),"qual5":set(),},
        "bcom":  {"qual0":set(),"qual1":set(),"qual2":set(),
                  "qual3":set(),"qual4":set(),"qual5":set(),},
      },
      "barnacle": {
        "total":set(),"lcom":set(),"gcom":set(),"bcom":set(),
      }
    }
    # barnacle_by_chrom[chromA][chromB][pos_keyA][pos_keyB] = list
    # chromosome "X" := "Xr" if align to "-" strand
    # pos_key = floor(pos / big_buffer)
    self.barnacle_by_chrom = dict()
    # barnacle_by_gene[geneA][geneB] = list
    # genes sorted by position
    self.barnacle_by_gene = dict()
    # barnacle contigs (used to determine primary breakpoints for
    #   multi-mapped predictions)
    self.multi_mapped_contigs = dict()
    # candidate_contig_map[full_candidate_id] = contig_id
    self.candidate_contig_map = dict()
    # defuse results to be aligned to contigs (qual2 or better)
    self.defuse_to_align = dict()
    self.defuse_to_realign = dict()
    # common gene-partners with tool-specific breakpoints
    self.breakpoint_partners = set()
    # breakpoints[chromA][chromB][pos_keyA][pos_keyB] = breakpoint_pair
    self.breakpoints = dict()
    self.outfiles = dict()
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
    if (hasattr(self, "outfiles")): #{
      for file_key in self.outfiles: #{
        self.outfiles[file_key].close()
      #} end for
    #} end if
  #} end def

  def RunComparison(self): #{
    LogMsg(self, "Comparing results...")
    start = time.time()
    # setup the output files
    self.SetupOutputFiles()
    # load Barnacle results
    self.LoadBarnacleResults()
    # process Defuse results
    self.ProcessDefuseResults()
    # count Barnacle events marked as common
    self.CountBarnacleCommons()
    if (self.options.align and 0 < len(self.defuse_to_align)): #{
      self.DoAlignmentAnalysis()
    #} end if
    # write count results
    self.WriteCounts()
    if (self.options.align): #{
      self.WriteAlignCounts()
    #} end if
    # write the breakpoint analysis
    self.WriteBreakpointAnalysis()
    LogMsg(self, "Time spent comparing results: %s" % TimeSpent(start))
  #} end def

  def SetupOutputFiles(self): #{
    LogMsg(self, "Preparing output files...")
    self.OpenOutputFile('barnacle', "only_barnacle.data",
      "unique Barnacle results")
    # only put qual3 or higher predictions in only Defuse file
    self.OpenOutputFile('defuse', "only_defuse.data", "unique Defuse results",
      header=DEFUSE_COMPARISON_HEADER)
    self.OpenOutputFile('counts', "comparison_counts.txt", "counts output",
      header="\n".join((COUNTS_HEADER1,COUNTS_HEADER2)))
    self.OpenOutputFile('defuse_seqs', "defuse_seqs.fa",
      "Defuse sequences output")
    self.OpenOutputFile('defuse_short_seqs', "defuse_short_seqs.fa",
      "Defuse short sequences output")
    self.OpenOutputFile('u_breakpoints', "unique_breakpoints.txt",
      "common gene-partners with tool-specific breakpoints")
    if (self.options.align): #{
      #self.OpenOutputFile('align', "defuse_to_contigs.psl",
      #  "Defuse-to-contig alignment results", mode="r")
      self.OpenOutputFile('align_results', "defuse_align_results.txt",
        "Processed results of Defuse-to-contig alignments",
        header=ALIGN_RESULTS_HEADER)
    #} end if
    if (self.options.realign): #{
      self.OpenOutputFile('unaligned_defuse_seqs', "unaligned_defuse_seqs.fa",
        "unaligned Defuse sequences output")
      #self.OpenOutputFile('realign', "defuse_to_unmerged.psl",
      #  "Defuse-to-unmerged contigs realignment results", mode="r")
      #self.OpenOutputFile('realign_results', "defuse_realign_results.txt",
      #  "Processed results of Defuse-to-unmerged contig realignments",
      #  header=ALIGN_RESULTS_HEADER)
    #} end if
    self.OpenOutputFile("breakpoints", "breakpoints.txt",
      "breakpoint analysis file")
  #} end def

  def OpenOutputFile(self, key, filename, desc, header=None, mode="w"): #{
    output_path = os.path.join(self.options.output_dir, filename)
    self.outfiles[key] = FileBoxCls(output_path, mode,
      "could not create %s file" % desc)
    if (None != header): #{
      self.outfiles[key].WriteLine(header)
    #} end if
  #} end def

  def LoadBarnacleResults(self): #{
    LogMsg(self, "Loading Barnacle results...")
    start = time.time()
    parser = CandidateGroupParserCls(self.options.barnacle_path)
    for event in parser: #{
      DebugMsg(self, "Loading event %i" % (event.id))
      for member in event.members: #{
        # track the Barnacle breakpoint
        #DebugMsg(self, "Getting breakpoint for member %s:\n"
        #  "  %s" % (member.IDString(), member.BreakpointsString()))
        member.id = member.IDString()
        matching_bp = self.GetBreakpoint(member)
        #DebugMsg(self, "  Adding Barnacle prediction to breakpoint:")
        matching_bp.AddBarnaclePrediction(member)
        #DebugMsg(self, "    %s" % matching_bp.ToString())
        ctg_id = member.contig_info.id
        if (0 == member.meta_fields['multi_mapped'] or
            ctg_id not in self.multi_mapped_contigs): #{
          matching_bp.primary = True
        #} end if
        if (0 < member.meta_fields['multi_mapped']): #{
          if (ctg_id not in self.multi_mapped_contigs): #{
            self.multi_mapped_contigs[ctg_id] = matching_bp
          #} end if
          matching_bp.multi_maps = True
          self.candidate_contig_map[member.IDString()] = ctg_id
        #} end if
      #} end for
      # presume event is not common
      event.lcommon = False
      event.gcommon = False
      event.bcommon = False
      # use the breakpoint from the first member
      breakpoint0 = event.members[0].breakpointA
      #breakpoint0.strand = event.members[0].align_info_A.Strand()
      breakpoint1 = event.members[0].breakpointB
      #breakpoint1.strand = event.members[0].align_info_B.Strand()
      (breakpoint0.genes, breakpoint1.genes) = event.GetGeneSets()
      sorted_breakpoints = SortBreakpoints(breakpoint0, breakpoint1)
      (event.breakpointA, event.breakpointB) = sorted_breakpoints
      #DebugMsg(self, "  %s\n  %s" %
      #  (event.breakpointA.ToString(), event.breakpointB.ToString()))
      # barnacle_by_chrom[chromA][chromB][pos_keyA][pos_keyB] = list
      # chromosome "X" := "Xr" if align to "-" strand
      chromA = event.breakpointA.chr
      chromB = event.breakpointB.chr
      #chromA = GetDirectedChromKey(event.breakpointA)
      #chromB = GetDirectedChromKey(event.breakpointB)
      # pos_key = floor(pos / big_buffer)
      posA = event.breakpointA.GetPositionKey(self.options.big_buffer)
      posB = event.breakpointB.GetPositionKey(self.options.big_buffer)
      #DebugMsg(self, "  chromA: %s, chromB: %s, posA: %i, posB: %i" %
      #  (chromA, chromB, posA, posB))
      SetupMultiDict(self.barnacle_by_chrom, [chromA, chromB, posA, posB],
        event)
      # barnacle_by_gene[geneA][geneB] = list
      # genes sorted by position
      for geneA in event.breakpointA.genes: #{
        for geneB in event.breakpointB.genes: #{
          #genesA = ",".join(sorted(event.breakpointA.genes))
          #genesB = ",".join(sorted(event.breakpointB.genes))
          SetupMultiDict(self.barnacle_by_gene, [geneA, geneB], event)
        #} end for
      #} end for
      #chromA_dict = SetupMultiDict(self.barnacle_by_chrom, chromA)
      #chromB_dict = SetupMultiDict(chromA_dict, chromB)
      #posA_dict = SetupMultiDict(chromB_dict, posA)
      #posB_list = SetupMultiDict(posA_dict, posB, bottom=True)
      #posB_list.append(event)
    #} end for
    if (self.log_info['debug']): #{
      partners = list()
      for genesA in self.barnacle_by_gene: #{
        for genesB in self.barnacle_by_gene[genesA]: #{
          partners.append(";".join([genesA, genesB]))
        #} end for
      #} end for
      LogMsg(self, "Barnacle gene partners:\n%s" % " ".join(sorted(partners)))
    #} end if
    LogMsg(self, "Time spent loading Barnacle results: %s" % TimeSpent(start))
  #} end def

  def GetBreakpoint(self, prediction): #{
    #DebugMsg(self, "  Getting breakpoint for prediction...")
    pair_tuple = SortBreakpoints(prediction.breakpointA,
      prediction.breakpointB)
    posA = pair_tuple[0].GetPositionKey(self.options.small_buffer)
    posB = pair_tuple[1].GetPositionKey(self.options.small_buffer)
    main_keys = [pair_tuple[0].chr, pair_tuple[1].chr, posA, posB]
    potential_bps = self.GetPotentialBreakpoints(main_keys)
    matching_bp = None
    #DebugMsg(self, "  Looking for existing breakpoint...")
    for breakpoint in potential_bps: #{
      #DebugMsg(self, "    checking breakpoint: %s" % breakpoint.ToString())
      if (BreakpointPairsAreEqual(pair_tuple, breakpoint.pair,
          buffer=breakpoint.buffer)): #{
        if (None == matching_bp): #{
          matching_bp = breakpoint
        else:
          LogMsg(self, "  Extra matching breakpoint found for %s\n"
            "  %s\n  %s" % (prediction.id, matching_bp.ToString(),
            breakpoint.ToString()))
        #} end if
      else:
        DebugMsg(self, "Breakpoints do not match:\n  %s\n  %s-%s" %
          (breakpoint.ToString(), pair_tuple[0].ToString(),
          pair_tuple[1].ToString()))
      #} end if
    #} end for
    if (None == matching_bp): #{
      #DebugMsg(self, "  Adding new breakpoint...")
      matching_bp = BreakpointPairCls(self.options.small_buffer)
      matching_bp.primary = False
      matching_bp.multi_maps = False
      try: #{
        SetupMultiDict(self.breakpoints, main_keys, matching_bp,
          list_end=False)
      except MultiDictError, e:
        clashing_bp = AccessMultiDict(self.breakpoints, main_keys)
        raise DefuseComparisonError("Clashing breakpoints with keys "
          "%s:\n    %s\n    %s" % (",".join(map(str,main_keys)),
          matching_bp.ToString(), clashing_bp.ToString()))
      #} end try
    #} end if
    return matching_bp
  #} end def

  def GetPotentialBreakpoints(self, main_keys): #{
    potential_breakpoints = list()
    # breakpoints[chromA][chromB][pos_keyA][pos_keyB] = breakpoint_pair
    for offsetA in [-1,0,1]: #{
      for offsetB in [-1,0,1]: #{
        keys = [main_keys[0], main_keys[1],
          main_keys[2]+offsetA, main_keys[3]+offsetB]
        if (CheckMultiDict(self.breakpoints, keys)): #{
          curr_breakpoint = AccessMultiDict(self.breakpoints, keys)
          potential_breakpoints.append(curr_breakpoint)
        #} end if
      #} end for
    #} end for
    return potential_breakpoints
  #} end def

  def ProcessDefuseResults(self): #{
    LogMsg(self, "Processing Defuse results...")
    start = time.time()
    defuse_file = FileBoxCls(self.options.defuse_path, "r",
      "could not read Defuse results file")
    field_cols = self.GetDefuseFieldColumns(defuse_file)
    for defuse_line in defuse_file: #{
      defuse_prediction = DefusePredictionCls(defuse_line, field_cols,
        log_info=self.log_info)
      self.ProcessDefusePrediction(defuse_prediction)
    #} end for
    LogMsg(self, "Time spent processing Defuse results: %s" % TimeSpent(start))
  #} end def

  def GetDefuseFieldColumns(self, defuse_file): #{
    LogMsg(self, "  Determining Defuse field order...")
    desired_fields = ["cluster_id", "splitr_sequence", "splitr_count",
      "adjacent", "altsplice", "deletion", "gene1", "gene2",
      "gene_chromosome1", "gene_chromosome2", "gene_location1",
      "gene_location2", "gene_name1", "gene_name2", "genomic_break_pos1",
      "genomic_break_pos2", "interchromosomal", "inversion", "library_name",
      "read_through", "probability",]
    field_cols = dict()
    heading_line = defuse_file.next()
    for (field_col, field_name) in enumerate(heading_line.split()): #{
      #DebugMsg(self, "  %i) %s" % (field_col, field_name))
      if (field_name in desired_fields): #{
        if (field_name in field_cols): #{
          raise DefuseComparisonError("non-unique field name: %s: %i, %i" %
            (field_name, field_col, field_cols[field_name]))
        #} end if
        field_cols[field_name] = field_col
      #} end if
    #} end for
    # check that desired fields were all found
    missing = list()
    for field_name in desired_fields: #{
      if (field_name not in field_cols): #{
        missing.append(field_name)
      #} end if
    #} end for
    if (0 < len(missing)): #{
      LogMsg(self, "WARNING missing Defuse fields: %s" % ", ".join(missing))
    #} end if
    return field_cols
  #} end def

  def ProcessDefusePrediction(self, defuse_prediction): #{
    DebugMsg(self, "  Processing Defuse prediction %s (%s/%s)..." % (
      defuse_prediction.id, defuse_prediction.breakpointA.gene_name,
      defuse_prediction.breakpointB.gene_name))
    # set the quality
    self.SetDefuseQuality(defuse_prediction)
    DebugMsg(self, "      Quality: %i" % defuse_prediction.quality)
    # update the total counts
    self.UpdateDefuseCounts("total", defuse_prediction.quality)
    self.UpdateDefusePartners("total", defuse_prediction)
    #bp_common = False
    if (self.LocationCommon(defuse_prediction)): #{
      self.UpdateDefuseCounts("lcom", defuse_prediction.quality)
      self.UpdateDefusePartners("lcom", defuse_prediction)
      for barnacle_pred in defuse_prediction.barnacle_lcom: #{
        barnacle_pred.lcommon = True
      #} end for
      if (self.PartnerCommon(defuse_prediction)): #{
        self.UpdateDefuseCounts("gcom", defuse_prediction.quality)
        self.UpdateDefusePartners("gcom", defuse_prediction)
      #} end if
      if (self.BreakpointCommon(defuse_prediction)): #{
        #bp_common = True
        self.UpdateDefuseCounts("bcom", defuse_prediction.quality)
        self.UpdateDefusePartners("bcom", defuse_prediction)
      elif (0 < len(defuse_prediction.barnacle_gcom)):
        partners = ";".join([defuse_prediction.breakpointA.gene_name,
          defuse_prediction.breakpointB.gene_name])
        self.breakpoint_partners.add("%s:defuse" % partners)
      #} end if
    elif (self.AltGeneCommon(defuse_prediction)): #{
      LogMsg(self, "Defuse prediction common by gene partners, but not by "
        "general location: %s: %s/%s" % (defuse_prediction.id,
        defuse_prediction.breakpointA.gene_name,
        defuse_prediction.breakpointB.gene_name))
      #common = True
      self.UpdateDefuseCounts("gcom", defuse_prediction.quality)
      self.UpdateDefusePartners("gcom", defuse_prediction)
      for barnacle_pred in defuse_prediction.barnacle_gcom: #{
        barnacle_pred.gcommon = True
      #} end for
    #} end if
    # track the Defuse breakpoint
    matching_bp = self.GetBreakpoint(defuse_prediction)
    matching_bp.AddOtherPrediction(defuse_prediction)
    # if the Defuse prediction is not common
    #if (not common and 2 <= defuse_prediction.quality): #{
    if (matching_bp.OtherOnly()): #{
      matching_bp.primary = True
      # if the Defuse prediction is not mitochondrial and not a readthrough
      if (2 <= defuse_prediction.quality): #{
        if (self.options.align): #{
          # store the prediction for when processing alignments
          # of Defuse sequences to assembled contigs
          if (-1 == defuse_prediction.breakpoint): #{
            LogMsg(self, "WARNING: could not get breakpoint coordinate for "
              "Defuse prediction %s" % defuse_prediction.id)
          else:
            self.defuse_to_align[defuse_prediction.id] = defuse_prediction
          #} end if
          # write Defuse sequence to file for aligning Defuse sequences
          # to assembled contigs
          self.outfiles['defuse_seqs'].WriteLine(">d%s\n%s" %
            (defuse_prediction.id, defuse_prediction.seq.replace("|","")))
          try: #{
            short_sequence = defuse_prediction.ShortSequence(
              self.DistFromBreakpoint())
          except DefusePredictionError, e:
            LogMsg(self, "Warning: %s" % e)
            short_sequence = defuse_prediction.seq
          #} end try
          self.outfiles['defuse_short_seqs'].WriteLine(">d%s\n%s" %
            (defuse_prediction.id, short_sequence))
        #} end if
        # if the Defuse prediction has a high quality
        if (3 <= defuse_prediction.quality): #{
          # write the Defuse prediction to the "only Defuse" output file
          DebugMsg(self, "    Recording good Defuse only prediction.")
          self.outfiles['defuse'].WriteLine(defuse_prediction.ToString())
        #} end if
      #} end if
    elif (matching_bp.Common() and not matching_bp.primary):
      for barnacle_id in matching_bp.barnacle_preds: #{
        contig_id = self.candidate_contig_map[barnacle_id]
        old_primary = self.multi_mapped_contigs[contig_id]
        old_primary.primary = False
        self.multi_mapped_contigs[contig_id] = matching_bp
      #} end for
    #} end if
  #} end def

  def SetDefuseQuality(self, defuse_pred): #{
    # qual0: total
    # qual1: not mitochondrial
    # qual2: not mitochondrial, not a readthrough-like event
    # qual3: not mitochondrial, not readthrough, and medium or better
    #         confidence (based on min value for 'probability' field)
    # qual4: not mitochondrial, not readthrough, and high confidence
    #         (based on min value for 'probability' field)
    # qual5: not mitochondrial, not readthrough, high confidence, and
    #         good coverage
    DebugMsg(self, "    Setting Defuse prediction quality")
    if ("M" in [defuse_pred.breakpointA.chr, defuse_pred.breakpointB.chr]): #{
      DebugMsg(self, "    Mitochondrial prediction!")
      return
    #} end if
    defuse_pred.quality += 1 # qual1
    if (defuse_pred.readthrough): #{
      DebugMsg(self, "    Readthrough-like prediction!")
      return
    #} end if
    defuse_pred.quality += 1 # qual2
    if (self.options.probability_threshold1 >= defuse_pred.prob): #{
      DebugMsg(self, "    Low confidence prediction (prob: %.2f)!" %
        defuse_pred.prob)
      return
    #} end if
    defuse_pred.quality += 1 # qual3
    if (self.options.probability_threshold2 >= defuse_pred.prob): #{
      DebugMsg(self, "    Medium confidence prediction (prob: %.2f)!" %
        defuse_pred.prob)
      return
    #} end if
    defuse_pred.quality += 1 # qual4
    if (self.options.coverage_threshold >= defuse_pred.coverage): #{
      DebugMsg(self, "    Low coverage prediction (coverage: %i)!" %
        defuse_pred.coverage)
      return
    #} end if
    defuse_pred.quality += 1 # qual5
    #if (defuse_pred.deletion): #{
    #  DebugMsg(self, "    Deletion prediction!")
    #  return
    ##} end if
    #defuse_pred.quality += 1 # qual5
  #} end def

  def UpdateDefuseCounts(self, count_type, pred_quality): #{
    for curr_qual in range(pred_quality+1): #{
      qual_key = "qual%i" % curr_qual
      self.counts['defuse'][count_type][qual_key] += 1
    #} end for
  #} end def

  def UpdateDefusePartners(self, count_type, pred): #{
    partners = ";".join(sorted([pred.breakpointA.gene_name,
      pred.breakpointB.gene_name]))
    for curr_qual in range(pred.quality+1): #{
      qual_key = "qual%i" % curr_qual
      self.partners['defuse'][count_type][qual_key].add(partners)
    #} end for
  #} end def

  def LocationCommon(self, defuse_pred): #{
    DebugMsg(self, "    Checking whether prediction is common by location")
    defuse_pred.barnacle_lcom = list()
    # barnacle_by_chrom[chromA][chromB][pos_keyA][pos_keyB] = list
    # pos_key = floor(pos / big_buffer)
    chromA = defuse_pred.breakpointA.chr
    chromB = defuse_pred.breakpointB.chr
    posA = defuse_pred.breakpointA.GetPositionKey(self.options.big_buffer)
    posB = defuse_pred.breakpointB.GetPositionKey(self.options.big_buffer)
    DebugMsg(self, "    chromA: %s, chromB: %s, posA: %i, posB: %i" %
      (chromA, chromB, posA, posB))
    if (not CheckMultiDict(self.barnacle_by_chrom, [chromA, chromB])): #{
      DebugMsg(self, "      Chromosome pair not common")
      return False
    #} end if
    for offsetA in [-1,0,1]: #{
      for offsetB in [-1,0,1]: #{
        keys = [chromA, chromB, posA+offsetA, posB+offsetB]
        DebugMsg(self, "    Keys: %s" % ",".join(map(str, keys)))
        if (not CheckMultiDict(self.barnacle_by_chrom, keys)): #{
          #DebugMsg(self, "      Keys not common.")
          continue
        #} end if
        barnacle_list = AccessMultiDict(self.barnacle_by_chrom, keys)
        if (0 == offsetA + offsetB): #{
          DebugMsg(self, "      Common predictions")
          defuse_pred.barnacle_lcom.extend(barnacle_list)
        else:
          for barnacle_pred in barnacle_list: #{
            lcommon = self.BreakpointsMatch(barnacle_pred, defuse_pred,
              self.options.big_buffer)
            if (lcommon): #{
              DebugMsg(self, "      Common predictions")
              defuse_pred.barnacle_lcom.append(barnacle_pred)
            #} end if
          #} end for
        #} end if
      #} end for
    #} end for
    return (0 < len(defuse_pred.barnacle_lcom))
  #} end def

  def BreakpointsMatch(self, barnacle_pred, defuse_pred, buffer): #{
    diffA = abs(barnacle_pred.breakpointA.coord -
      defuse_pred.breakpointA.coord)
    diffB = abs(barnacle_pred.breakpointB.coord -
      defuse_pred.breakpointB.coord)
    DebugMsg(self, "      Comparing breakpoint coordinates (buffer:%s)\n"
      "      Barnacle: %s; %s\n      Defuse: %s; %s\n      DiffA:%i DiffB:%i" %
      (buffer, barnacle_pred.breakpointA.ToString(),
       barnacle_pred.breakpointB.ToString(),
       defuse_pred.breakpointA.ToString(),
       defuse_pred.breakpointB.ToString(), diffA, diffB))
    if (buffer < diffA or buffer < diffB): #{
      return False
    #} end if
    return True
  #} end def

  def PartnerCommon(self, defuse_pred): #{
    DebugMsg(self, "    Checking whether prediction is common by "
      "gene-partners")
    defuse_pred.barnacle_gcom = list()
    dgeneA = defuse_pred.breakpointA.gene_name
    dgeneB = defuse_pred.breakpointB.gene_name
    DebugMsg(self, "    Defuse: %s; %s" % (dgeneA, dgeneB))
    for barnacle_pred in defuse_pred.barnacle_lcom: #{
      DebugMsg(self, "    Barnacle: %s; %s" % (
        ",".join(sorted(barnacle_pred.breakpointA.genes)),
        ",".join(sorted(barnacle_pred.breakpointB.genes))))
      if (dgeneA in barnacle_pred.breakpointA.genes and
          dgeneB in barnacle_pred.breakpointB.genes): #{
        DebugMsg(self, "      Common predictions")
        defuse_pred.barnacle_gcom.append(barnacle_pred)
        barnacle_pred.gcommon = True
      #} end if
    #} end for
    return (0 < len(defuse_pred.barnacle_gcom))
  #} end def

  def BreakpointCommon(self, defuse_pred): #{
    DebugMsg(self, "    Checking whether prediction is common by breakpoint")
    defuse_pred.barnacle_bcom = list()
    for barnacle_pred in defuse_pred.barnacle_lcom: #{
      bcommon = self.BreakpointsMatch(barnacle_pred, defuse_pred,
        self.options.small_buffer)
      if (bcommon): #{
        DebugMsg(self, "      Common predictions")
        defuse_pred.barnacle_bcom.append(barnacle_pred)
        barnacle_pred.bcommon = True
      #} end if
    #} end for
    return (0 < len(defuse_pred.barnacle_bcom))
  #} end def

  def AltGeneCommon(self, defuse_pred): #{
    DebugMsg(self, "    Checking whether prediction is common by "
      "gene-partners but not location")
    defuse_pred.barnacle_gcom = list()
    keys = [defuse_pred.breakpointA.gene_name,
      defuse_pred.breakpointB.gene_name]
    DebugMsg(self, "    Keys: %s" % ",".join(map(str, keys)))
    if (CheckMultiDict(self.barnacle_by_gene, keys)): #{
      DebugMsg(self, "    Common partners!")
      barnacle_list = AccessMultiDict(self.barnacle_by_gene, keys)
      defuse_pred.barnacle_gcom.extend(barnacle_list)
      for barnacle_pred in barnacle_list: #{
        if (defuse_pred.breakpointA.chr != barnacle_pred.breakpointA.chr or
            defuse_pred.breakpointB.chr != barnacle_pred.breakpointB.chr): #{
          LogMsg(self, "WARNING predictions common by genes, but not by "
            "chromosomes!\nBarnacle: %s; %s\nDefuse: %s; %s" %
            (barnacle_pred.breakpointA.ToString(),
             barnacle_pred.breakpointB.ToString(),
             defuse_pred.breakpointA.FullString(),
             defuse_pred.breakpointB.FullString()))
        #} end if
      #} end for
    #} end if
    return (0 < len(defuse_pred.barnacle_gcom))
  #} end def

  def CountBarnacleCommons(self): #{
    LogMsg(self, "Counting common Barnacle predictions...")
    start = time.time()
    # use barnacle_by_chrom
    for barnacle_pred in TraverseMultiDict(self.barnacle_by_chrom): #{
      partners = ";".join([",".join(sorted(barnacle_pred.breakpointA.genes)),
        ",".join(sorted(barnacle_pred.breakpointB.genes))])
      DebugMsg(self, "  G%i: %s lcom:%s gcom:%s bcom:%s" %
        (barnacle_pred.id, partners, barnacle_pred.lcommon,
        barnacle_pred.gcommon, barnacle_pred.bcommon))
      self.UpdateBarnacleCounts("total")
      gene_partners = (barnacle_pred.breakpointA.genes,
        barnacle_pred.breakpointB.genes)
      self.UpdateBarnaclePartners("total", gene_partners)
      #DebugMsg(self, "Barnacle gene partners:\n%s" %
      #  " ".join((";".join([",".join(sorted(partner[0])),
      #                      ",".join(sorted(partner[1]))]) for
      #  partner in self.partners['barnacle']['total'])))
      common = False
      if (barnacle_pred.lcommon): #{
        common = True
        self.UpdateBarnacleCounts("lcom")
        self.UpdateBarnaclePartners("lcom", gene_partners)
        #DebugMsg(self, "Barnacle lcommon gene partners:\n%s" %
        #  " ".join((";".join([",".join(sorted(partner[0])),
        #                      ",".join(sorted(partner[1]))]) for
        #  partner in self.partners['barnacle']['lcom'])))
      #} end if
      if (barnacle_pred.gcommon): #{
        common = True
        self.UpdateBarnacleCounts("gcom")
        self.UpdateBarnaclePartners("gcom", gene_partners)
        #DebugMsg(self, "Barnacle gcommon gene partners:\n%s" %
        #  " ".join((";".join([",".join(sorted(partner[0])),
        #                      ",".join(sorted(partner[1]))]) for
        #  partner in self.partners['barnacle']['gcom'])))
      #} end if
      if (barnacle_pred.bcommon): #{
        common = True
        self.UpdateBarnacleCounts("bcom")
        self.UpdateBarnaclePartners("bcom", gene_partners)
        #DebugMsg(self, "Barnacle bcommon gene partners:\n%s" %
        #  " ".join((";".join([",".join(sorted(partner[0])),
        #                      ",".join(sorted(partner[1]))]) for
        #  partner in self.partners['barnacle']['bcom'])))
      elif (barnacle_pred.gcommon):
        self.breakpoint_partners.add("%s:barnacle" % partners)
      #} end if
      # if the Barnacle prediction is not common
      if (not common): #{
        DebugMsg(self, "    Recording Barnacle only prediction.")
        # write the Barnacle prediction to the "only Barnacle" output file
        self.outfiles['barnacle'].Write(barnacle_pred.FullDataString())
      #} end if
    #} end for
    LogMsg(self, "Time spent counting common Barnacle predictions: %s" %
      TimeSpent(start))
  #} end def

  def UpdateBarnacleCounts(self, count_type): #{
    self.counts['barnacle'][count_type] += 1
  #} end def

  def UpdateBarnaclePartners(self, count_type, new_partners): #{
    DebugMsg(self, "Updating Barnacle %s gene-partners" % count_type)
    for old_partners in self.partners['barnacle'][count_type]: #{
      commonA = old_partners[0].intersection(new_partners[0])
      commonB = old_partners[1].intersection(new_partners[1])
      if (0 < len(commonA) and 0 < len(commonB)): #{
        if (commonA == new_partners[0] and
            commonB == new_partners[1]): #{
          DebugMsg(self, "New partners already present")
          return
        #} end if
        DebugMsg(self, "Updating old partners %s with new partners %s" % (
          ";".join([",".join(sorted(old_partners[0])),
                    ",".join(sorted(old_partners[1]))]),
          ";".join([",".join(sorted(new_partners[0])),
                    ",".join(sorted(new_partners[1]))])))
        new_partners[0].update(old_partners[0])
        new_partners[1].update(old_partners[1])
        self.partners['barnacle'][count_type].remove(old_partners)
        break
      #} end if
    #} end for
    frozen_partners = (frozenset(new_partners[0]), frozenset(new_partners[1]))
    self.partners['barnacle'][count_type].add(frozen_partners)
  #} end def

  def DoAlignmentAnalysis(self): #{
    # align Defuse sequences to Contigs
    self.AlignDefuseSequences()
    # process Defuse-to-Contig alignments
    self.ProcessAlignments()
    if (self.options.realign): #{
      # select the Defuse predictions to realign
      self.SelectDefusePredictionsToRealign()
      if (0 < len(self.defuse_to_realign)): #{
        for kvalue in range(self.options.low_k, self.options.high_k+1): #{
          self.DoRealignmentAnalysis(kvalue=kvalue, contig_type="contigs")
          if (self.options.include_bubbles): #{
            self.DoRealignmentAnalysis(kvalue=kvalue, contig_type="bubbles")
          #} end if
          if (self.options.include_indels): #{
            self.DoRealignmentAnalysis(kvalue=kvalue, contig_type="indel")
          #} end if
        #} end for
      #} end if
    #} end if
  #} end def

  def DoRealignmentAnalysis(self, kvalue, contig_type): #{
    key = "%i-%s" % (kvalue, contig_type)
    # align Defuse sequences to Contigs
    self.AlignDefuseSequences(key)
    # process Defuse-to-Contig alignments
    self.ProcessAlignments(key)
    # collate the information and reset for next k-value
    self.FinishProcessingK(key)
  #} end def

  def AlignDefuseSequences(self, key=None): #{
    start = time.time()
    # get the appropriate target, query, and output paths
    if (None == key): #{
      LogMsg(self, "Aligning Defuse sequences to assembled contigs...")
      target_path = self.options.contigs_path
      query_file = self.outfiles['defuse_short_seqs']
      output_path = self.options.align_path
    else:
      LogMsg(self, "Realigning unaligned Defuse sequences to unmerged "
        "k%s..." % key)
      target_path = self.options.realign_target_paths[key]
      query_file = self.outfiles['unaligned_defuse_seqs']
      output_path = self.options.realign_output_paths[key]
    #} end if
    # close the Defuse sequences query file
    query_file.Close()
    # use "-minMatch=1", "stepSize=5", "minScore=4" to allow short matches
    # use "repMatch=2253" to mimic webBLAT results
    command = [GetCommand(self, "blat"), "-minMatch=1", "-stepSize=5",
      "-repMatch=2253", "-minScore=4", target_path, query_file.path,
      output_path]
    DebugMsg(self, " ".join(command))
    return_code = RunCommandFromList(command, dpt=self.options.dpt)
    if (0 != return_code and not self.log_info['debug']): #{
      LogMsg(self, " ".join(command))
    #} end if
    if (0 > return_code): #{
      raise DefuseComparisonError("Alignment command was terminated by "
        "signal %i" % return_code)
    elif (0 < return_code):
      raise DefuseComparisonError("Error running alignment command: %i" %
        return_code)
    #} end if
    LogMsg(self, "Time spent aligning Defuse sequences: %s" % TimeSpent(start))
  #} end def

  def ProcessAlignments(self, key=None): #{
    LogMsg(self, "Processing Defuse-to-contig alignments...")
    start = time.time()
    # get the appropriate list of defuse predictions
    if (None == key): #{
      preds_to_process = self.defuse_to_align
      aligns_path = self.options.align_path
    else:
      preds_to_process = self.defuse_to_realign
      aligns_path = self.options.realign_output_paths[key]
    #} end if
    filters = { 'identity': 80.0 }
    try:
      aligns = ParseAlignmentFile(aligns_path, filters, log_info=self.log_info)
    except NoAlignmentsError, e:
      LogMsg(self, "  Found no alignments")
      return
    #} end if
    for align in map(FixAlign, aligns): #{
      ExtremeDebugMsg(self, "Processing alignment %s->%s: %i-%i, M:%i" %
        (align.query, align.target, align.qstart, align.qend, align.match))
      # skip the 'd' at the beginning of the query name
      cluster_id = align.query[1:]
      if (cluster_id not in preds_to_process): #{
        ExtremeDebugMsg(self, "  skipping alignment for prediction %s" %
          cluster_id)
        continue
      #} end if
      defuse_pred = preds_to_process[cluster_id]
      if (3 > align.qstart and 2 > (align.query_len - align.qend) and
          3 > align.qbaseinsert and 3 > align.tbaseinsert): #{
        DebugMsg(self, "Found full Defuse alignment!\n  %s" % align.psl_str)
        defuse_pred.full_align = True
        #defuse_pred.UpdateBestAligns(align)
        defuse_pred.best_aligns['full'].Update(align)
        continue
      #} end if
      breakpoint = self.DistFromBreakpoint()
      overlapA = (breakpoint - align.qstart) + 1
      overlapB = (align.qend - breakpoint) + 1
      ExtremeDebugMsg(self, "  Cluster %s, breakpoint: %i overlap: %i/%i" %
        (defuse_pred.ShortString(), breakpoint, overlapA, overlapB))
      if (self.options.align_threshold2 <= min(overlapA, overlapB)): #{
        overlaps = self.GetTrueOverlaps(breakpoint, align)
        tleft = min(align.tstart, align.tend)
        tright = max(align.tstart, align.tend)
        if (3 > tleft and 2 > (align.target_len - tright) and
            3 > align.tbaseinsert and
            self.options.align_threshold2 <= min(overlaps)): #{
          DebugMsg(self, "Found full contig alignment!\n  %s" %
            align.psl_str)
          defuse_pred.full_align = True
          #defuse_pred.UpdateBestAligns(align)
          defuse_pred.best_aligns['tfull'].Update(align)
        elif (self.options.align_threshold1 <= min(overlaps)):
          defuse_pred.part_align = True
          defuse_pred.best_aligns['part'].Update(align)
          #if (not defuse_pred.full_align): #{
          #  defuse_pred.UpdateBestAligns(align)
          #} end if
        #} end if
      #} end if
    #} end for
    LogMsg(self, "Time spent processing alignments: %s" % TimeSpent(start))
  #} end def

  def GetTrueOverlaps(self, breakpoint, align): #{
    ExtremeDebugMsg(self, "    Getting true overlap values")
    # overlaps = [total_before_bp, total_after_bp,
    #   continuous_before_bp, continuous_after_bp]
    overlaps = [0,0,0,0]
    (dprev,cprev) = (None,None)
    continuous_after = True
    for (defuse_block,ctg_block) in zip(align.query_blocks,
        align.target_blocks): #{
      #(dblock_start, dblock_end) = sorted(map(int, defuse_block))
      #dsize = (dblock_end-dblock_start)+1
      dblock = AlignBlockCls(defuse_block)
      cblock = AlignBlockCls(ctg_block)
      if (dblock.span != cblock.span): #{
        raise DefuseComparisonError("defuse sequence and contig alignment "
          "block sizes do not match: %i, %i" % (dblock.span, cblock.span))
      #} end if
      if (None != dprev and None != cprev): #{
        if (dblock.start < breakpoint and
            (CONTINUOUS_THRESHOLD > dblock.Gap(dprev) or
             CONTINUOUS_THRESHOLD > cblock.Gap(cprev))): #{
          overlaps[2] = 0
        #} end if
        if (breakpoint < dblock.start and
            (CONTINUOUS_THRESHOLD > dblock.Gap(dprev) or
             CONTINUOUS_THRESHOLD > cblock.Gap(cprev))): #{
          continuous_after = False
        #} end if
      #} end if
      if (dblock.end < breakpoint): #{
        overlaps[0] += dblock.span
        overlaps[2] += dblock.span
      elif (breakpoint < dblock.start):
        overlaps[1] += dblock.span
        if (continuous_after): #{
          overlaps[3] += dblock.span
        #} end if
      else:
        before = (breakpoint-dblock.start)+1
        overlaps[0] += before
        overlaps[2] += before
        after = (dblock.end-breakpoint)+1
        overlaps[1] += after
        overlaps[3] += after
      #} end if
      ExtremeDebugMsg(self, "Defuse block: %s, Contig block: %s, Overlaps: "
        "%s (%s)" % (dblock, cblock, "/".join(map(str,overlaps[:2])),
        "/".join(map(str, overlaps[2:]))))
      (dprev, cprev) = (dblock, cblock)
    #} end for
    return overlaps
  #} end def

  def SelectDefusePredictionsToRealign(self): #{
    for pred in self.defuse_to_align.itervalues(): #{
      # if the defuse prediction sequence did not align to any merged contigs
      if (not pred.full_align and not pred.part_align): #{
        pred.realigned = True
        self.defuse_to_realign[pred.id] = pred
        # write Defuse sequence to file for realigning Defuse sequences
        # to pre-merging assembled contigs
        try: #{
          short_sequence = pred.ShortSequence(self.DistFromBreakpoint())
        except DefusePredictionError, e:
          LogMsg(self, "Warning: %s" % e)
          short_sequence = defuse_prediction.seq
        #} end try
        self.outfiles['unaligned_defuse_seqs'].WriteLine(">d%s\n%s" %
          (pred.id, short_sequence))
      #} end if
    #} end for
  #} end def

  def FinishProcessingK(self, key): #{
    for pred in self.defuse_to_realign.itervalues(): #{
      if (pred.full_align): #{
        pred.any_full_realign = True
      #} end if
      if (pred.part_align): #{
        pred.any_part_realign = True
      #} end if
      # set the appropriate values for the current k-value and clear the
      # general results to prepare for the next k-value
      pred.full_realign[key] = pred.full_align
      pred.full_align = False
      pred.part_realign[key] = pred.part_align
      pred.part_align = False
      pred.best_realigns[key] = dict()
      for align_type in ['full', 'tfull', 'part']: #{
        align_list = pred.best_aligns[align_type].aligns
        pred.best_realigns[key][align_type] = align_list
        pred.best_aligns[align_type].Reset()
      #} end for
    #} end for
  #} end def

  def WriteCounts(self): #{
    #LogMsg(self, " ".join(self.partners['defuse']['total']['qual0']))
    LogMsg(self, "%s\n%s" % (COUNTS_HEADER1, COUNTS_HEADER2))
    for qual in range(6): #{
      qual_key = "qual%i" % qual
      data_list = ["defuse_%i" % qual]
      for count_type in ["total","lcom","gcom","bcom"]: #{
        data_list.append(self.counts['defuse'][count_type][qual_key])
      #} end for
      for count_type in ["total","lcom","gcom","bcom"]: #{
        data_list.append(len(self.partners['defuse'][count_type][qual_key]))
      #} end for
      data_str = "\t".join(map(str, data_list))
      self.outfiles['counts'].WriteLine(data_str)
      LogMsg(self, data_str)
    #} end for
    data_list = ["barnacle"]
    for count_type in ["total","lcom","gcom","bcom"]: #{
      data_list.append(self.counts['barnacle'][count_type])
    #} end for
    for count_type in ["total","lcom","gcom","bcom"]: #{
      data_list.append(len(self.partners['barnacle'][count_type]))
    #} end for
    data_str = "\t".join(map(str, data_list))
    self.outfiles['counts'].WriteLine(data_str)
    LogMsg(self, "%s\n%s" % (self.options.lib, data_str))
    # write out common gene-partners with tool-specific breakpoints
    for partners in sorted(self.breakpoint_partners): #{
      self.outfiles['u_breakpoints'].WriteLine(partners)
    #} end for
  #} end def

  def WriteAlignCounts(self): #{
    total = len(self.defuse_to_align)
    full_aligns = part_aligns = no_aligns = 0
    realign_total = len(self.defuse_to_realign)
    full_realigns = part_realigns = no_realigns = 0
    for defuse_pred in self.defuse_to_align.itervalues(): #{
      if (defuse_pred.realigned): #{
        no_aligns += 1
        if (defuse_pred.any_full_realign): #{
          full_realigns += 1
        elif (defuse_pred.any_part_realign):
          part_realigns += 1
        else:
          no_realigns += 1
        #} end if
      else:
        if (defuse_pred.full_align): #{
          full_aligns += 1
        elif (defuse_pred.part_align):
          part_aligns += 1
        else:
          no_aligns += 1
        #} end if
      #} end if
      self.outfiles['align_results'].WriteLine(defuse_pred.AlignString())
    #} end for
    data_str = ("Attempted to align %s Defuse sequences\n  Full: %i, "
      "Partial: %i, None: %i" % (total, full_aligns, part_aligns, no_aligns))
    if (self.options.realign): #{
      data_str += ("\nAttempted to realign %s Defuse sequences\n  Full: %i, "
        "Partial: %i, None: %i" % (realign_total, full_realigns, part_realigns,
        no_realigns))
    #} end if
    self.outfiles['counts'].WriteLine(data_str)
    LogMsg(self, data_str)
  #} end def

  def WriteBreakpointAnalysis(self): #{
    DebugMsg(self, "Writing breakpoint analysis...")
    num_common = 0
    num_common_mm = 0
    num_barnacle = 0
    num_barnacle_mm = 0
    num_defuse = 0
    breakpoint_iter = TraverseMultiDict(self.breakpoints, sort=True)
    for (index, breakpoint) in enumerate(breakpoint_iter): #{
      if (breakpoint.Common()): #{
        DebugMsg(self, "Breakpoint %i is common" % index)
        if (breakpoint.primary): #{
          num_common += 1
          if (breakpoint.multi_maps): #{
            num_common_mm += 1
          #} end if
        else:
          LogMsg(self, "Breakpoint %i is common, but not primary!" % index)
        #} end if
      elif (breakpoint.BarnacleOnly()):
        DebugMsg(self, "Breakpoint %i is Barnacle-specific" % index)
        if (breakpoint.primary): #{
          DebugMsg(self, "  and primary!")
          num_barnacle += 1
          if (breakpoint.multi_maps): #{
            num_barnacle_mm += 1
          #} end if
        #} end if
      elif (breakpoint.OtherOnly()):
        DebugMsg(self, "Breakpoint %i is Defuse-specific" % index)
        # only count high quality Defuse specific breakpoints
        if (None == breakpoint.other_quality or
            3 <= breakpoint.other_quality): #{
          DebugMsg(self, "  and good quality!")
          num_defuse += 1
        #} end if
      else:
        LogMsg(self, "MYSTERY BREAKPOINT! %s" % breakpoint.ToString())
      #} end if
      DebugMsg(self, "  %s" % breakpoint.ToString())
      self.outfiles['breakpoints'].WriteLine("%i %s" % (index,
        breakpoint.ToString()))
    #} end for
    #self.outfiles['breakpoints'].WriteLine("\nCommon breakpoints: %i" %
    #  num_common)
    #self.outfiles['breakpoints'].WriteLine("Barnacle specific breakpoints: "
    #  "%i" % num_barnacle)
    #self.outfiles['breakpoints'].WriteLine("Defuse specific breakpoints "
    #  "with qual > 2: %i" % num_defuse)
    self.outfiles['breakpoints'].WriteLine("\n".join(["",
      "Common breakpoints: %i (%i multi-map)" % (num_common, num_common_mm),
      "Barnacle specific breakpoints: %i (%i multi-map)" %
        (num_barnacle, num_barnacle_mm),
      "Defuse specific breakpoints with qual > 2: %i" % num_defuse,
    ]))
  #} end def

  def DistFromBreakpoint(self): #{
    return int(self.options.align_threshold1 * 1.5)
  #} end def
#} end class

def GetDirectedChromKey(breakpoint): #{
  # chromosome "X" := "Xr" if align to "-" strand
  chrom_key = breakpoint.chr
  if ("-" == breakpoint.strand): #{
    chrom_key += "r"
  #} end if
  return chrom_key
#} end def

#### EXCEPTION CLASSES ####
class DefuseComparisonError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("") #TODO
  args = [ "LIB", "BARNACLE_FILE", "DEFUSE_FILE", "CONTIGS_FILE" ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-B", "--big-buffer",
                    type="int", metavar="N",
                    help="Consider events similar if their breakpoint "
                      "coordinates are within Nbp of each other. "
                      "[default:%default]")
  parser.add_option("-b", "--small-buffer",
                    type="int", metavar="N",
                    help="Consider events identical if their breakpoint "
                      "coordinates are within Nbp of each other. "
                      "[default:%default]")
  parser.add_option("-p", "--probability-threshold1",
                    type="float", metavar="F",
                    help="Consider Defuse events to be \"medium confidence\" "
                      "if the probability field has a value greater than F. "
                      "[default:%default]")
  parser.add_option("-P", "--probability-threshold2",
                    type="float", metavar="F",
                    help="Consider Defuse events to be \"high confidence\" if "
                      "the probability field has a value greater than F. "
                      "[default:%default]")
  parser.add_option("-c", "--coverage-threshold",
                    type="int", metavar="N",
                    help="Consider Defuse events to have \"good coverage\" if "
                      "the splitr_count field has a value greater than N. "
                      "[default:%default]")
  parser.add_option("-a", "--align",
                    action="store_true",
                    help="Align Defuse prediction sequences with unique "
                      "breakpoints to assembled contigs.")
  parser.add_option("-A", "--align-threshold1",
                    type="int", metavar="N",
                    help="Consider a Defuse-to-contig alignment \"partial\" "
                      "if less than the full Defuse sequence is aligned, but "
                      "the alignment overlaps the breakpoint by at least Nbp "
                      "on either side. [default:%default]")
  parser.add_option("--align-threshold2",
                    type="int", metavar="N",
                    help="Consider a Defuse-to-contig alignment \"partial\" "
                      "if less than the full Defuse sequence is aligned, but "
                      "the full contig is aligned and it overlaps the "
                      "breakpoint by at least Nbp on either side. "
                      "[default:%default]")
  parser.add_option("--realign-unaligned",
                    action="store_true", dest="realign",
                    help="Realign Defuse prediction sequences with unique "
                      "breakpoints that do not align to the merged contig "
                      "assembly to a pre-merging assembly.")
  parser.add_option("--realign-base",
                    dest="realign_dir", metavar="DIR",
                    help="The base directory holding pre-merging assembled "
                      "contig sequence files to use during realignment.")
  parser.add_option("--low-k",
                    type="int",
                    help="The minimum k-value pre-merging assembly to use "
                      "when realigning. [default: %default]")
  parser.add_option("--high-k",
                    type="int",
                    help="The maximum k-value pre-merging assembly to use "
                      "when realigning. [default: %default]")
  parser.add_option("--include-bubbles",
                    action="store_true",
                    help="When realigning unique Defuse prediction sequences "
                      "to pre-merging assembled contigs, also align to "
                      "\"-bubbles.fa\" file. [default]")
  parser.add_option("--exclude-bubbles",
                    action="store_false", dest="include_bubbles",
                    help="When realigning unique Defuse prediction sequences "
                      "to pre-merging assembled contigs, do not align to "
                      "\"-bubbles.fa\" file.")
  parser.add_option("--include-indels",
                    action="store_true",
                    help="When realigning unique Defuse prediction sequences "
                      "to pre-merging assembled contigs, also align to "
                      "\"-indels.fa\" file. [default]")
  parser.add_option("--exclude-indels",
                    action="store_false", dest="include_indels",
                    help="When realigning unique Defuse prediction sequences "
                      "to pre-merging assembled contigs, do not align to "
                      "\"-indels.fa\" file.")
  parser.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force comparison to take place, even if the output "
                         "directory already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.set_defaults(big_buffer=100000,
                      small_buffer=50,
                      probability_threshold1=0.75,
                      probability_threshold2=0.9,
                      coverage_threshold=3,
                      align=False,
                      align_threshold1=23,
                      align_threshold2=19,
                      realign=False,
                      low_k=26,
                      high_k=50,
                      include_bubbles=True,
                      include_indels=True,
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
  if (options.align_threshold1 < options.align_threshold2): #{
    ErrMsg("Alignment threshold 1 must be larger than alignment threshold 2")
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle results",  path_errors)
  CheckFilePath(options.defuse_path,   "Defuse results",    path_errors)
  CheckFilePath(options.contigs_path,  "assembled contigs", path_errors)
  # get and check the output path
  options.output_dir = GetOutDir(os.path.dirname(options.barnacle_path),
    "defuse_comparison")
  if (options.realign): #{
    # need to do original alignment if realigning
    options.align = True
    if (options.low_k > options.high_k): #{
      ErrMsg("minimum k-value (%i) must not be greater than maximum k-value "
        "(%i)" % (options.low_k, options.high_k))
      opts_good = False
    #} end if
    if (None == options.realign_dir): #{
      ErrMsg("you must provide the path to a directory holding pre-merging "
        "assembly contig sequence files (--realign-base) when using the "
        "\"--realign\" option.")
      opts_good = False
    else:
      CheckDirPath(options.realign_dir, "pre-merging assembly base",
        path_errors)
    #} end if
    realign_out_dir = os.path.join(options.output_dir, "unaligned_realign")
    options.realign_target_paths = dict()
    options.realign_output_paths = dict()
    for kvalue in range(options.low_k, options.high_k+1): #{
      SetupRealignPath(options, realign_out_dir, kvalue, "contigs",
        path_errors)
      if (options.include_bubbles): #{
        SetupRealignPath(options, realign_out_dir, kvalue, "bubbles",
          path_errors)
      #} end if
      if (options.include_indels): #{
        SetupRealignPath(options, realign_out_dir, kvalue, "indel",
          path_errors)
      #} end if
    #} end for
  #} end if
  if (options.align): #{
    options.align_path = os.path.join(options.output_dir,
      "defuse_to_contigs.psl")
  #} end if
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "comparison output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "defuse_comparison", options.output_dir)
  #} end if
  if (opts_good and 0 == len(path_errors) and options.realign): #{
    CheckDirPath(realign_out_dir, "realignment output", path_errors,
      create=True, replace=True)
  #} end if
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # the paths are good if there are no path errors and no conflicting options
  return (opts_good and 0 == len(path_errors))
#} end def

def SetupRealignPath(options, realign_out_dir, kvalue, contig_type,
    path_errors): #{
  key = "%i-%s" % (kvalue, contig_type)
  options.realign_target_paths[key] = os.path.join(options.realign_dir,
    "k%i" % kvalue, "%s-%s.fa" % (options.lib, contig_type))
  CheckFilePath(options.realign_target_paths[key],
    "k%i %s assembly" % (kvalue, contig_type), path_errors)
  options.realign_output_paths[key] = os.path.join(realign_out_dir,
    "defuse_to_k%i-%s.psl" % (kvalue, contig_type))
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.lib           = args[0]
    options.barnacle_path = EnsureAbsPath(args[1])
    options.defuse_path   = EnsureAbsPath(args[2])
    options.contigs_path  = EnsureAbsPath(args[3])
    if (CheckPaths(options)): #{
      try: #{
        compare_boss = DefuseComparisonCls(options)
        WriteCommand(compare_boss, sys.argv)
        compare_boss.RunComparison()
      except (MyError), e:
        ErrMsg("ERROR while comparing results:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle data file for that library (BARNACLE_FILE); the path to a "
      "Defuse results file for that library (DEFUSE_FILE); and the path to "
      "a fasta-formatted assembled contig sequences file (CONTIGS_FILE).")
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
