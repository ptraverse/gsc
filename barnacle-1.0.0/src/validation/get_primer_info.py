#! /usr/bin/env python
"""
get_primer_info.py

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
  ReverseComplement)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls, CheckNewFilePath)
from parsers.fasta import FastaFileCls
from parsers.candidate_group_parser import CandidateGroupParserCls
from alignment_processing.alignment_functions import (ParseAlignmentFile,
  FixAlign)
from primer_info import PrimerInfoCls, PrimerSeqCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"
MIN_CTG_LEN = 1000

class PrimerDesignerCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    # self.primer_infos[ctg_id] = primer_info_list
    self.primer_infos = dict()
    # self.contigs[ctg_id] = sequence
    self.contigs = dict()
    # self.transcripts[gene_id] = sequence
    self.transcripts = dict()
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def GetPrimerInfo(self): #{
    LogMsg(self, "Designing primers...")
    start = time.time()
    # load the data needed to design primers
    self.LoadData()
    # create the output files
    self.CreateOutputFiles()
    # write the info
    self.WriteInfo()
    # write the sequences
    self.WriteSequences()
    LogMsg(self, "Time spent designing primers: %s" % TimeSpent(start))
  #} end def

  def LoadData(self): #{
    LogMsg(self, "Loading data...")
    start = time.time()
    # load the prediction data
    LogMsg(self, "Loading Barnacle predictions...")
    group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    for group in group_parser: #{
      type_is_valid = False
      for valid_type in ["fusion", "ptd", "itd"]: #{
        if (valid_type in group.event_type): #{
          type_is_valid = True
          break
        #} end if
      #} end for
      if (not type_is_valid): #{
        raise PrimerDesignerError("unsupported event type: %s" %
          group.event_type)
      #} end if
      member = self.SelectMember(group)
      DebugMsg(self, "Selected member:\n%s" % member.DataString())
      self.PreProcessMember(member, group.event_type)
    #} end for
    # load the alignments for the chosen members
    for aligns_path in self.options.aligns_paths: #{
      self.LoadAlignments(aligns_path)
    #} end for
    for primer_info_list in self.primer_infos.itervalues(): #{
      for primer_info in primer_info_list: #{
        DebugMsg(self, "PRIMER_INFO: %s\nALIGNS_0: %s (%i)\n"
          "ALIGNS_1: %s (%i)" % (primer_info.ToString(),
          primer_info.c2t_aligns[0], len(primer_info.c2t_aligns_extra[0]),
          primer_info.c2t_aligns[1], len(primer_info.c2t_aligns_extra[1])))
      #} end for
    #} end for
    # load the contig sequences
    self.LoadContigSequences()
    # load the required transcript sequences
    self.LoadTranscriptSequences()
    LogMsg(self, "Time spent loading data: %s" % TimeSpent(start))
  #} end def

  def SelectMember(self, group): #{
    DebugMsg(self, "Selecting member for group %i with %i members" %
      (group.id, len(group.members)))
    if (1 == len(group.members)): #{
      return group.members[0]
    #} end if
    max_ctg_len = 0
    max_short_align = 0
    for member in group.members: #{
      max_ctg_len = max(max_ctg_len, member.contig_info.length)
      member.short_align = min(member.align_info_A.ContigSpan(),
        member.align_info_B.ContigSpan())
      max_short_align = max(max_short_align, member.short_align)
    #} end for
    selected_member = None
    DebugMsg(self, "MAX contig length: %i, MAX shorter alignment: %i" %
      (max_ctg_len, max_short_align))
    for member in group.members: #{
      ctg_len_score = float(member.contig_info.length) / float(max_ctg_len)
      short_align_score = float(member.short_align) / float(max_short_align)
      member.primer_score = ctg_len_score + short_align_score
      DebugMsg(self, "Member %s: contig length score = %.2f, shorter "
        "alignment score = %.2f, total score = %.2f" % (member.IDString(),
        ctg_len_score, short_align_score, member.primer_score))
      if (None == selected_member or
          selected_member.primer_score < member.primer_score): #{
        selected_member = member
      #} end if
    #} end for
    DebugMsg(self, "SELECTED %s: total score = %.2f" %
      (selected_member.IDString(), selected_member.primer_score))
    return selected_member
  #} end def

  def PreProcessMember(self, member, event_type): #{
    primer_info = PrimerInfoCls(member, event_type, self.log_info)
    DebugMsg(self, "PRIMER_INFO: %s" % primer_info.ToString())
    if (primer_info.ctg_id not in self.contigs):
      self.contigs[primer_info.ctg_id] = None
    #} end if
    if (primer_info.ctg_id not in self.primer_infos):
      self.primer_infos[primer_info.ctg_id] = list()
    #} end if
    self.primer_infos[primer_info.ctg_id].append(primer_info)
  #} end def

  def LoadAlignments(self, aligns_path): #{
    filters = { 'count': 100, 'identity': 90.0 }
    aligns = ParseAlignmentFile(aligns_path, filters, log_info=self.log_info)
    LogMsg(self, "Loading alignment data...")
    start = time.time()
    for align in map(FixAlign, aligns): #{
      ctg_id = align.query
      if (ctg_id not in self.primer_infos): #{
        continue
      #} end if
      #DebugMsg(self, "Loading alignment of %s" % ctg_id)
      for primer_info in self.primer_infos[ctg_id]: #{
        primer_info.UpdateC2TAligns(align)
        #DebugMsg(self, "PRIMER_INFO: %s\nALIGNS_0: %s (%i)\n"
        #  "ALIGNS_1: %s (%i)" % (primer_info.ToString(),
        #  primer_info.c2t_aligns[0], len(primer_info.c2t_aligns_extra[0]),
        #  primer_info.c2t_aligns[1], len(primer_info.c2t_aligns_extra[1])))
      #} end for
    #} end for
    LogMsg(self, "Time spent loading alignment data: %s" % TimeSpent(start))
  #} end def

  def LoadContigSequences(self): #{
    LogMsg(self, "Loading contig sequences...")
    start = time.time()
    ctg_seq_file = FastaFileCls(self.options.ctg_seq_path,
      fail_msg="cannot read contig sequences path")
    for seq_object in ctg_seq_file: #{
      if (seq_object.id not in self.contigs): #{
        continue
      #} end if
      if (self.contigs[seq_object.id] != None and
          self.contigs[seq_object.id] != seq_object.sequence): #{
        raise PrimerDesignerError("conflicting sequences for %s contig:\n"
          "%s\n%s" % (seq_object.id, seq_object.sequence,
          self.contigs[seq_object.id]))
      #} end if
      self.contigs[seq_object.id] = seq_object.sequence
      #DebugMsg(self, seq_object.Output())
    #} end for
    LogMsg(self, "Time spent loading contig sequences: %s" % TimeSpent(start))
  #} end def

  def LoadTranscriptSequences(self): #{
    LogMsg(self, "Loading transcript sequences...")
    start = time.time()
    for primer_info_list in self.primer_infos.itervalues(): #{
      for primer_info in primer_info_list: #{
        primer_info.ChooseAlignments()
        self.transcripts[primer_info.c2t_aligns[0].transcript] = None
        if ("itd" not in primer_info.event_type): #{
          self.transcripts[primer_info.c2t_aligns[1].transcript] = None
        #} end if
      #} end for
    #} end for
    transcript_file = FastaFileCls(self.options.tran_seq_path,
      fail_msg="cannot read transcript sequences path")
    for seq_object in transcript_file: #{
      if (seq_object.id not in self.transcripts): #{
        continue
      #} end if
      if (self.transcripts[seq_object.id] != None and
          self.transcripts[seq_object.id] != seq_object.sequence): #{
        raise PrimerDesignerError("conflicting sequences for %s "
          "transcript:\n%s\n%s" % (seq_object.id, seq_object.sequence,
          self.transcripts[seq_object.id]))
      #} end if
      self.transcripts[seq_object.id] = seq_object.sequence
      #DebugMsg(self, seq_object.Output())
    #} end for
    LogMsg(self, "Time spent loading transcript sequences: %s" %
      TimeSpent(start))
  #} end def

  def CreateOutputFiles(self): #{
    DebugMsg(self, "Creating primer info file: %s" % self.options.info_path)
    self.info_file = FileBoxCls(self.options.info_path, "w",
      "cannot create primer info output file")
    DebugMsg(self, "Creating sequence file: %s" % self.options.seqs_path)
    self.seqs_file = FileBoxCls(self.options.seqs_path, "w",
      "cannot create primer design sequences file")
  #} end def

  def WriteInfo(self): #{
    LogMsg(self, "Writing primer info...")
    start = time.time()
    for primer_info_list in self.primer_infos.itervalues(): #{
      for primer_info in primer_info_list: #{
        self.OutputPrimerInfo(primer_info)
      #} end for
    #} end for
    LogMsg(self, "Time spent writing primer info: %s" % TimeSpent(start))
  #} end def

  def OutputPrimerInfo(self, primer_info): #{
    if ("fusion" in primer_info.event_type): #{
      data_list = self.GetFusionPrimerInfo(primer_info)
    elif ("ptd" in primer_info.event_type): #{
      data_list = self.GetPTDPrimerInfo(primer_info)
    elif ("itd" in primer_info.event_type): #{
      data_list = self.GetITDPrimerInfo(primer_info)
    else:
      raise PrimerDesignerError("unrecognized event type: %s" %
        primer_info.event_type)
    #} end if
    data_str = "\t".join(data_list)
    DebugMsg(self, "PRIMER_INFO: %s\nALIGNS_0: %s (%i)\n"
      "ALIGNS_1: %s (%i)\n%s" % (primer_info.ToString(),
      primer_info.c2t_aligns[0], len(primer_info.c2t_aligns_extra[0]),
      primer_info.c2t_aligns[1], len(primer_info.c2t_aligns_extra[1]),
      data_str))
    self.info_file.WriteLine(data_str)
  #} end def

  def GetFusionPrimerInfo(self, primer_info): #{
    # [event_id, contig_id, contig_left, contig_right, geneA_id,
    #  geneA_left, geneA_right, geneB_id, geneB_left, geneB_right]
    # event_id
    event_id = primer_info.EventID()
    # contig_id
    contig = PrimerSeqCls(primer_info.ctg_id)
    # contig_left
    contig.left = min(primer_info.c2t_aligns[0].ctg_end,
      primer_info.c2t_aligns[1].ctg_start)
    # contig_right
    contig.right = max(primer_info.c2t_aligns[0].ctg_end,
      primer_info.c2t_aligns[1].ctg_start)
    # contig_overlap
    ctg_overlap = primer_info.C2TAlignOverlap()
    # ensure that geneA is the 5' gene, if possible
    if ("+" == primer_info.c2t_aligns[0].strand): #{
      alignA = primer_info.c2t_aligns[0]
      alignB = primer_info.c2t_aligns[1]
      geneA = self.GeneInfoFromAlign0(alignA, ctg_overlap)
      geneB = self.GeneInfoFromAlign1(alignB, ctg_overlap)
    else:
      alignA = primer_info.c2t_aligns[1]
      alignB = primer_info.c2t_aligns[0]
      geneA = self.GeneInfoFromAlign1(alignA, ctg_overlap)
      geneB = self.GeneInfoFromAlign0(alignB, ctg_overlap)
    #} end if
    if (MIN_CTG_LEN > len(self.contigs[contig.id])): #{
      DebugMsg(self, "Contig %s too short (%i): attempting to construct "
        "synthetic chimera" % (contig.id, len(self.contigs[contig.id])))
      chimera = self.ConstructFusionSeq(primer_info, contig, ctg_overlap)
      if (chimera.id in self.transcripts and
          self.transcripts[chimera.id] != chimera.seq): #{
        raise PrimerDesignerError("conflicting chimeric fusion sequences: %s" %
          chimera.id)
      #} end if
      if (len(chimera.seq) > len(self.contigs[contig.id])): #{
        DebugMsg(self, "Using synthetic chimera sequence")
        self.transcripts[chimera.id] = chimera.seq
        contig = chimera
      #} end if
    #} end if
    return map(str, [event_id, contig.id, contig.left, contig.right,
      geneA.id, geneA.left, geneA.right, geneB.id, geneB.left, geneB.right])
  #} end def

  def GeneInfoFromAlign0(self, align, overlap): #{
    gene = PrimerSeqCls(align.transcript)
    # use the end of the alignment
    gene.left = gene.right = align.transcript_end
    if (0 < overlap): #{
      if ("+" == align.strand): #{
        gene.left -= overlap
      else:
        gene.right += overlap
      #} end if
    #} end if
    return gene
  #} end def

  def GeneInfoFromAlign1(self, align, overlap): #{
    gene = PrimerSeqCls(align.transcript)
    # use the start of the alignment
    gene.left = gene.right = align.transcript_start
    if (0 < overlap): #{
      if ("+" == align.strand): #{
        gene.right += overlap
      else:
        gene.left -= overlap
      #} end if
    #} end if
    return gene
  #} end def

  def ConstructFusionSeq(self, primer_info, contig, ctg_overlap): #{
    # chimera_id
    chimera = PrimerSeqCls("G%i_%s_%s/%s" % (primer_info.group_id,
      primer_info.ctg_id, primer_info.c2t_aligns[0].transcript,
      primer_info.c2t_aligns[1].transcript))
    left_seq  = self.RefSeqFromAlign0(primer_info.c2t_aligns[0])
    right_seq = self.RefSeqFromAlign1(primer_info.c2t_aligns[1])
    extra = ""
    if (0 > ctg_overlap): #{
      extra = self.contigs[primer_info.ctg_id][contig.left:contig.right-1]
    else:
      right_seq = right_seq[ctg_overlap:]
    #} end if
    chimera.seq = (left_seq + extra + right_seq)
    # chimera_left
    chimera.left = min(len(left_seq), len(left_seq)+ctg_overlap)
    # chimera_right
    chimera.right = max(len(left_seq), len(left_seq)+ctg_overlap)
    DebugMsg(self, "Left: %s\nExtra: %s\nRight: %s\nOverlap: %i, "
      "WT_Length0: %i, WT_Length1: %i, Extra_len: %i, CT_Length: %i, "
      "Calc_Length: %i" % (left_seq, extra, right_seq, ctg_overlap,
      len(left_seq), len(right_seq)+max(0,ctg_overlap), len(extra),
      len(chimera.seq), len(left_seq)+len(right_seq)+len(extra)))
    return chimera
  #} end def

  def RefSeqFromAlign0(self, align): #{
    start_extra = ""
    if (1 < align.ctg_start): #{
      if (1 == align.transcript_start or
          align.transcript_len == align.transcript_start): #{
        start_extra = self.contigs[align.ctg_id][:align.ctg_start-1]
      #} end if
    #} end if
    DebugMsg(self, "Ctg Len: %i, Ctg Start: %i, Trn Len: %i, Trn Start: %i\n"
      "Extra: %s" % (align.ctg_len, align.ctg_start, align.transcript_len,
      align.transcript_start, start_extra))
    # use the end of the alignment
    ref_seq = self.transcripts[align.transcript]
    coord = align.transcript_end
    if ("+" == align.strand): #{
      return start_extra + ref_seq[:coord]
    else:
      return start_extra + ReverseComplement(ref_seq[coord-1:])
    #} end if
  #} end def

  def RefSeqFromAlign1(self, align): #{
    end_extra = ""
    if (align.ctg_len > align.ctg_end): #{
      if (1 == align.transcript_end or
          align.transcript_len == align.transcript_end): #{
        end_extra = self.contigs[align.ctg_id][align.ctg_end:]
      #} end if
    #} end if
    DebugMsg(self, "Ctg Len: %i, Ctg End: %i, Trn Len: %i, Trn End: %i\n"
      "Extra: %s" % (align.ctg_len, align.ctg_end, align.transcript_len,
      align.transcript_end, end_extra))
    # use the start of the alignment
    ref_seq = self.transcripts[align.transcript]
    coord = align.transcript_start
    if ("+" == align.strand): #{
      return ref_seq[coord-1:] + end_extra
    else:
      return ReverseComplement(ref_seq[:coord]) + end_extra
    #} end if
  #} end def

  def GetPTDPrimerInfo(self, primer_info): #{
    # [event_id, chimera.id, chimera.left, chimera.right, gene.id,
    #  gene.left, gene.right]
    # event_id
    event_id = primer_info.EventID()
    # contig_id
    contig = PrimerSeqCls(primer_info.ctg_id)
    # contig_left
    contig.left = min(primer_info.c2t_aligns[0].ctg_end,
      primer_info.c2t_aligns[1].ctg_start)
    # contig_right
    contig.right = max(primer_info.c2t_aligns[0].ctg_end,
      primer_info.c2t_aligns[1].ctg_start)
    # gene_id
    gene = PrimerSeqCls(primer_info.c2t_aligns[0].transcript)
    # gene.left
    gene.left = min(primer_info.c2t_aligns[0].transcript_end,
      primer_info.c2t_aligns[1].transcript_start)
    # gene.right
    gene.right = max(primer_info.c2t_aligns[0].transcript_end,
      primer_info.c2t_aligns[1].transcript_start)
    # dup_length
    dup_length = (gene.right - gene.left) + 1
    # construct inferred chimera sequence
    ctg_overlap = primer_info.C2TAlignOverlap()
    extra = ""
    if (0 > ctg_overlap): #{
      extra = self.contigs[primer_info.ctg_id][contig_left+1:contig_right]
    #} end if
    if ("-" == primer_info.c2t_aligns[0].strand): #{
      extra = ReverseComplement(extra)
    #} end if
    # chimera_id
    chimera = PrimerSeqCls("G%i_%s_%s" % (primer_info.group_id,
      primer_info.ctg_id, gene.id))
    chimera.seq = (self.transcripts[gene.id][:gene.right] + extra +
      self.transcripts[gene.id][gene.left-1:])
    if (chimera.id in self.transcripts and
        self.transcripts[chimera.id] != chimera.seq): #{
      raise PrimerDesignerError("conflicting chimeric PTD sequences: %s" %
        chimera.id)
    #} end if
    self.transcripts[chimera.id] = chimera.seq
    # chimera_left
    chimera.left = gene.left
    # chimera_right
    chimera.right = gene.right + len(extra) + dup_length
    DebugMsg(self, "Overlap: %i, Extra: \"%s\", WT_Length: %i, Dup_Length: "
      "%i, CT_Length: %i, Calc_Length: %i" % (ctg_overlap, extra,
      len(self.transcripts[gene.id]), dup_length, len(chimera.seq),
      len(self.transcripts[gene.id])+len(extra)+dup_length))
    return map(str, [event_id, chimera.id, chimera.left, chimera.right,
      gene.id, gene.left, gene.right])
  #} end def

  def GetITDPrimerInfo(self, primer_info): #{
    # [event_id, chimera.id, chimera.left, chimera.right, gene.id,
    #  gene.left, gene.right]
    DebugMsg(self, "Getting Primer Info for %s %s" %
      (primer_info.EventID(), primer_info.ctg_id))
    # event_id
    event_id = primer_info.EventID()
    # contig_id
    contig = PrimerSeqCls(primer_info.ctg_id)
    # contig_left
    contig.left = min(primer_info.c2g_aligns[1][0],
      primer_info.c2g_aligns[2][0])
    # contig_right
    contig.right = max(primer_info.c2g_aligns[1][1],
      primer_info.c2g_aligns[2][1])
    # alignment to transcript
    align = primer_info.c2t_aligns[0]
    # ensure that there are some blocks
    if (None == align.query_blocks): #{
      align.query_blocks =[[align.ctg_start,align.ctg_end]]
    #} end if
    if (None == align.target_blocks): #{
      align.target_blocks =[
        [align.transcript_start,align.transcript_end]]
    #} end if
    #align = self.AdjustITDAlign(primer_info)
    #DebugMsg(self, "Adjusted: %s" % align.ToString())
    #cdup_coords = primer_info.c2g_aligns[1]
    cdup_coords = self.GetContigDuplicationCoords(primer_info.c2g_aligns,
      align)
    DebugMsg(self, "CDup: %i-%i" % (cdup_coords[0], cdup_coords[1]))
    # gene_id
    gene = PrimerSeqCls(align.transcript)
    num_blocks = len(align.query_blocks)
    for index in range(num_blocks): #{
      DebugMsg(self, "Block: %i-%i::%i-%i" % (align.query_blocks[index][0],
        align.query_blocks[index][1], align.target_blocks[index][0],
        align.target_blocks[index][1]))
      if (align.query_blocks[index][0] <= cdup_coords[0] and
          cdup_coords[1] <= align.query_blocks[index][1]): #{
        DebugMsg(self, "Found block %i" % index)
        break
      #} end if
    #} end for
    if (num_blocks <= index): #{
      raise PrimerDesignerError("could not find block containing duplication!")
    #} end if
    start_offset = (cdup_coords[0] - align.query_blocks[index][0])
    if ("-" == align.strand): #{
      start_offset = -start_offset
    #} end if
    tdup_start = align.target_blocks[index][0] + start_offset
    #for index in range(index, num_blocks): #{
    #  if (align.query_blocks[index][1] > cdup_coords[1]): #{
    #    break
    #  #} end if
    #} end for
    end_offset = (cdup_coords[1] - align.query_blocks[index][0])
    if ("-" == align.strand): #{
      end_offset = -end_offset
    #} end if
    tdup_end = align.target_blocks[index][0] + end_offset
    DebugMsg(self, "In Block: %i-%i::%i-%i, offsets:%i,%i, tdup:%i-%i" % (
      align.query_blocks[index][0], align.query_blocks[index][1],
      align.target_blocks[index][0], align.target_blocks[index][1],
      start_offset, end_offset, tdup_start, tdup_end))
    # gene.left
    gene.left = min(tdup_start, tdup_end)
    # gene.right
    gene.right = max(tdup_start, tdup_end)
    cdup_seq = self.contigs[contig.id][cdup_coords[0]-1:cdup_coords[1]]
    tdup_seq = self.transcripts[gene.id][gene.left-1:gene.right]
    if ("-" == align.strand): #{
      tdup_seq = ReverseComplement(tdup_seq)
    #} end if
    DebugMsg(self, "Align: %s\nDup-coords: %s-%s;%s-%s\n  %s\n  %s" % (
      align.ToString(), cdup_coords[0], cdup_coords[1], tdup_start,
      tdup_end, cdup_seq, tdup_seq))
    if (cdup_seq.upper() != tdup_seq.upper()): #{
      diffs = 0
      for i in range(len(cdup_seq)): #{
        if (cdup_seq[i].upper() != tdup_seq[i].upper()): #{
          diffs += 1
        #} end if
      #} end for
      if (1 < diffs): #{
        raise PrimerDesignerError("Contig duplicated-sequence is not the "
          "same as transcript duplicated-sequence!")
      else:
        LogMsg(self, "single-base mismatch from ref in duplicated seq")
      #} end if
    #} end if
    return map(str, [event_id, contig.id, contig.left, contig.right,
      gene.id, gene.left, gene.right])
  #} end def

  def AdjustITDAlign(self, primer_info): #{
    align = primer_info.c2t_aligns[0]
    DebugMsg(self, "Adjusting alignment: %s\n  %s" %
      (primer_info.EventID(), align.ToString()))
    old_gap_coords = primer_info.c2g_aligns[2]
    gap_len = abs(old_gap_coords[0] - old_gap_coords[1]) + 1
    if (not primer_info.internal_gap): #{
      DebugMsg(self, "  gap is not internal!")
      return align
    #} end if
    cseq = self.contigs[primer_info.ctg_id]
    tseq = self.transcripts[align.transcript]
    if (None == align.query_blocks): #{
      align.query_blocks =[[align.ctg_start,align.ctg_end]]
    #} end if
    if (None == align.target_blocks): #{
      align.target_blocks =[[align.transcript_start,align.transcript_end]]
    #} end if
    if ("-" == align.strand): #{
      # adjust the blocks to refer to the reverse complement
      for block in align.target_blocks: #{
        block[0] = align.transcript_len - block[0] + 1
        block[1] = align.transcript_len - block[1] + 1
      #} end if
    #} end if
    # if the gap was moved to the start of the contig
    if (3 < align.ctg_start and align.ctg_start < old_gap_coords[0]): #{
      DebugMsg(self, "cseq: %s\ntseq: %s" % (cseq, tseq))
      i = align.ctg_start-1
      seq1 = cseq[i:i+gap_len]
      k = cseq.find(seq1)
      DebugMsg(self, "seq1: %s, gap:%i i:%i j:%i k:%i l:%i" % (seq1, gap_len,
        i, i+gap_len-1, k, k+gap_len-1))
      DebugMsg(self, "Comparing: %s--%s" % (cseq[i+gap_len], cseq[k+gap_len]))
      while (cseq[i+gap_len] == cseq[k+gap_len] and k+gap_len < i): #{
        gap_len += 1
        seq1 = cseq[i:i+gap_len]
        DebugMsg(self, "seq1: %s, gap:%i i:%i j:%i k:%i l:%i" % (seq1, gap_len,
          i, i+gap_len-1, k, k+gap_len-1))
        DebugMsg(self, "Comparing: %s--%s" %
          (cseq[i+gap_len], cseq[k+gap_len]))
      #} end while
      # split the first block
      if (k < i): #{
        old_qblock = align.query_blocks[0]
        new_qblock1 = [k+1,k+gap_len]
        new_qblock2 = [i+gap_len+1,old_qblock[1]]
        old_tblock = align.target_blocks[0]
        new_tblock1 = [old_tblock[0],old_tblock[0]+gap_len-1]
        new_tblock2 = [old_tblock[0]+gap_len,old_tblock[1]]
        align.query_blocks[:1] = [new_qblock1,new_qblock2]
        align.target_blocks[:1] = [new_tblock1,new_tblock2]
      #} end if
      extra = 0
      DebugMsg(self, "seq1: %s, extra:%i k':%i l:%i" % (seq1, extra,
        k-extra, k+gap_len-1))
      kt = align.transcript_start-1
      DebugMsg(self, "Comparing: %s--%s" % (cseq[k-extra-1], tseq[kt-extra-1]))
      while (extra < k and
          cseq[k-extra-1].upper() == tseq[kt-extra-1].upper()): #{
        extra += 1
        seq1 = cseq[k-extra:k+gap_len]
        DebugMsg(self, "seq1: %s, extra:%i k':%i l:%i" % (seq1, extra,
          k-extra, k+gap_len-1))
        DebugMsg(self, "Comparing: %s--%s" %
          (cseq[k-extra-1], tseq[kt-extra-1]))
      #} end while
      # extend the first block
      align.query_blocks[0][0] -= extra
      align.target_blocks[0][0] -= extra
      # shift the gap to the right
      shift = 0
      seq1 = cseq[i:i+gap_len].lower() + cseq[i+gap_len].upper()
      DebugMsg(self, "seq1: %s, shift:%i gap:%i-%i" % (seq1, shift,
        i+shift, i+gap_len+shift-1))
      DebugMsg(self, "Comparing: %s--%s" %
        (cseq[i+shift], cseq[i+gap_len+shift]))
      while (cseq[i+shift] == cseq[i+gap_len+shift]): #{
        shift += 1
        seq1 = (cseq[i+shift:i+gap_len+shift].lower() +
          cseq[i+gap_len+shift].upper())
        DebugMsg(self, "seq1: %s, shift:%i gap:%i-%i" % (seq1, shift,
          i+shift, i+gap_len+shift-1))
        DebugMsg(self, "Comparing: %s--%s" %
          (cseq[i+shift], cseq[i+gap_len+shift]))
      #} end while
      align.query_blocks[0][1] += shift
      align.query_blocks[1][0] += shift
      align.target_blocks[0][1] += shift
      align.target_blocks[1][0] += shift
    #} end if
    # if the gap was moved to the start of the contig
    if (3 < (align.ctg_len - align.ctg_end) and
        old_gap_coords[1] < align.ctg_end): #{
      DebugMsg(self, "cseq: %s\ntseq: %s" % (cseq, tseq))
      raise PrimerDesignerError("not implemented!")
    #} end if
    if ("-" == align.strand): #{
      # restore the blocks from the reverse complement
      for block in align.target_blocks: #{
        block[0] = align.transcript_len - block[0] + 1
        block[1] = align.transcript_len - block[1] + 1
      #} end if
    #} end if
    # update the contig start and end
    align.ctg_start = align.query_blocks[0][0]
    align.ctg_end   = align.query_blocks[-1][1]
    # update the transcript start and end
    align.transcript_start = align.target_blocks[0][0]
    align.transcript_end   = align.target_blocks[-1][1]
    return align
  #} end def

  def GetContigDuplicationCoords(self, c2g_aligns, c2t_align): #{
    # use the contig duplication coords overlapped by the c2t alignment
    for block in c2t_align.query_blocks: #{
      if (block[0] <= c2g_aligns[1][0] and c2g_aligns[1][1] <= block[1]): #{
        return c2g_aligns[1]
      #} end if
    #} end for
    return c2g_aligns[2]
  #} end def

  def WriteSequences(self): #{
    LogMsg(self, "Writing sequences...")
    start = time.time()
    for ctg_id in self.contigs: #{
      self.OutputSequence(ctg_id, self.contigs[ctg_id])
    #} end for
    for gene_id in self.transcripts: #{
      self.OutputSequence(gene_id, self.transcripts[gene_id])
    #} end for
    LogMsg(self, "Time spent writing sequences: %s" % TimeSpent(start))
  #} end def

  def OutputSequence(self, seq_id, sequence): #{
    self.seqs_file.WriteLine(">%s" % seq_id)
    self.seqs_file.WriteLine(sequence)
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class PrimerDesignerError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Use Primer3 and gfPcr to create qPCR and sequence "
    "validation primers for the Barnacle events in the given file.")
  args = [ "LIB", "BARNACLE_FILE", "ALIGNMENTS_FILES", "CONTIG_SEQS_FILE",
    "TRANSCRIPT_SEQS_FILE" ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(force=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors)
  for aligns_path in options.aligns_paths: #{
    CheckFilePath(aligns_path, "contig-to-transcriptome alignments",
      path_errors)
  #} end for
  CheckFilePath(options.ctg_seq_path, "contig sequences", path_errors)
  CheckFilePath(options.tran_seq_path, "transcript sequences", path_errors)
  # get and check the output path
  options.output_dir = os.path.join(os.path.dirname(options.barnacle_path),
    "validation")
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=True)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "validation", options.output_dir)
    options.info_path = os.path.join(options.output_dir,
      "%s.primer.info" % os.path.basename(options.barnacle_path))
    if (not options.force): #{
      CheckNewFilePath(options.info_path, "primer info", path_errors)
    #} end if
    options.seqs_path = os.path.join(options.output_dir, "%s.seqs.fa" %
      os.path.basename(options.barnacle_path))
    if (not options.force): #{
      CheckNewFilePath(options.seqs_path, "sequences", path_errors)
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
    options.aligns_paths  = [EnsureAbsPath(path) for
      path in args[2].split(",")]
    options.ctg_seq_path  = EnsureAbsPath(args[3])
    options.tran_seq_path = EnsureAbsPath(args[4])
    if (CheckPaths(options)): #{
      try: #{
        primer_designer = PrimerDesignerCls(options)
        WriteCommand(primer_designer, sys.argv)
        primer_designer.GetPrimerInfo()
      except (MyError), e:
        ErrMsg("ERROR while designing primers:\n  %s" % e) #TODO
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); the path to a "
      "Barnacle data file for that library (BARNACLE_FILE); a comma-delimited "
      "list of paths to psl-formatted candidate contig to transcript sequence "
      "alignment files (ALIGNMENTS_FILES); the path to a fasta-formatted candidate "
      "contig sequence file (CONTIG_SEQS_FILE); and the path to a fasta-formatted "
      "transcript sequence file (TRANSCRIPT_SEQS_FILE).")
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
