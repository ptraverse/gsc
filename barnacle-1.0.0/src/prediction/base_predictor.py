#! /usr/bin/env python
"""
base_predictor.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.general import SetupMainClass
from utils.messages import LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import GetFilePath, FileBoxCls
from parsers.fasta import FastaFileCls
from common.candidate_group import CandidateGroupCls
from common.grouped_candidate import GroupedCandidateCls
from parsers.candidate_group_parser import CandidateGroupParserCls

class BasePredictorCls:
  def __init__(self, key, ext, description, options, log_info=None):
    SetupMainClass(self, options, log_info=log_info)
    self.key = key
    self.ext = ext
    self.description = description
    self.write_pretty = True
    self.write_special = True
    if (options.realign): #{
      self.write_pretty = False
      self.write_special = False
    #} end if
    self.outfiles = dict()
    self.AddOutFile("main", self.ext, self.description, options.realign)
    #output_path = GetFilePath(self.output_dir, options.lib, self.ext)
    #self.output_file = FileBoxCls(output_path, "w",
    #  "cannot create %s output file" % self.description)
    if (options.pretty): #{
      self.AddOutFile("pretty", "%s.pretty" % self.ext, "pretty %s" %
        self.description)
      #pretty_path = GetFilePath(options.output_dir, options.lib,
      #  "%s.pretty" % self.ext)
      #self.pretty_file = FileBoxCls(pretty_path, "w",
      #  "cannot create pretty %s output file" % self.description)
    #else:
    #  self.pretty_file = None
    #} end if
    self.num_predictions = 0
    self.transcript_sequences = dict()
    # only certain event types require the sequence for the realignment test
    #self.store_seq = False
  # end def

  def AddOutFile(self, key, ext, desc, use_realign_dir=False): #{
    if (use_realign_dir): #{
      output_dir = self.options.realign_dir
    else:
      output_dir = self.options.output_dir
    #} end if
    output_path = GetFilePath(output_dir, self.options.lib, ext)
    self.outfiles[key] = FileBoxCls(output_path, "w",
      "cannot create %s output file" % desc)
  #} end def

  def __del__(self): #{
    #self.output_file.Close()
    for file in self.outfiles.itervalues(): #{
      file.close()
    #} end for
  #} end def

  def ProcessGroup(self, group, good_members=None): #{
    DebugMsg(self, "Testing whether group %i is a %s:" %
      (group.id, self.description.upper()))
    ExtremeDebugMsg(self, "%s" % group.DataString())
    # check whether group satisfies event-type criteria
    event = None
    if (self.TestGroup(group)): #{
      event = CandidateGroupCls(group.DataString())
    else:
      DebugMsg(self, "Group is not a %s" % self.description)
      return False
    #} end if
    for (index, member) in enumerate(group.members): #{
      DebugMsg(self, "  Testing whether candidate contig %s is a %s:" %
        (member.IDString(), self.description))
      ExtremeDebugMsg(self, "  %s" % member.DataString())
      # check whether member satisfies event-type criteria
      if (self.TestMember(member)): #{
        #ExtremeDebugMsg(self, "  AFTER TEST: %s" % member.DataString())
        event.members.append(GroupedCandidateCls(member.DataString()))
        if (None != good_members): #{
          good_members.append(index)
        #} end if
      else:
        DebugMsg(self, "  Member is not a %s" % self.description)
        continue
      #} end if
    #} end for
    #if (None != event): #{
    #  ExtremeDebugMsg(self, "AFTER TEST:\n%s" % event.FullDataString())
    #} end if
    event.UpdateMemberCount()
    # if no candidate contigs satisfy the event-type criteria
    if (0 == event.num_members): #{
      DebugMsg(self, "Group has no %s members" % self.description)
      return False
    #} end if
    DebugMsg(self, "Predicting %s event from group!" % self.description)
    event.event_type = self.key
    self.OutputEvent(event)
    return True
  #} end def

  def OutputEvent(self, event): #{
    self.num_predictions += 1
    # write the event to the output file(s)
    #self.output_file.Write(event.FullDataString())
    self.outfiles["main"].Write(event.FullDataString())
    #if (None != self.pretty_file and self.write_pretty): #{
    if ("pretty" in self.outfiles and self.write_pretty): #{
      #self.pretty_file.Write(event.FullDataString(readable=True))
      self.outfiles["pretty"].Write(event.FullDataString(readable=True))
    #} end if
  #} end def

  def LoadTranscriptSequences(self, realigned_contigs): #{
    transcripts_of_interest = set()
    for realigned_contig in realigned_contigs.itervalues(): #{
      if (None == realigned_contig.best_align): #{
        DebugMsg(self, "Contig %s has no best realignment" %
          realigned_contig.id)
      else:
        transcripts_of_interest.add(realigned_contig.best_align.target)
      #} end if
      for region_pair in realigned_contig.region_pairs.itervalues(): #{
        transcripts_of_interest.update(region_pair.gene_sets['A'])
        transcripts_of_interest.update(region_pair.gene_sets['B'])
      #} end for
      for full_dup_targets in realigned_contig.full_dup_aligns.values(): #{
        #ExtremeDebugMsg(self, "  Adding full dup target(s): %s" %
        #  ",".join(full_dup_targets))
        transcripts_of_interest.update(full_dup_targets)
      #} end for
    #} end for
    DebugMsg(self, "Will attempt to load %i transcript sequences" %
      len(transcripts_of_interest))
    t_seqs_file = FastaFileCls(self.options.tran_seq_path,
      fail_msg="cannot open transcript sequences file")
    num_seqs = num_loaded = 0
    for seq in t_seqs_file: #{
      if (seq.id in transcripts_of_interest): #{
        self.transcript_sequences[seq.id] = seq.sequence.lower()
        num_loaded += 1
      #} end if
      num_seqs += 1
    #} end for
    DebugMsg(self, "Processed %i sequences, loaded %i" %
      (num_seqs, num_loaded))
  #} end def

  #def ReprocessPredictions(self, contig_alignments, contig_sequences): #{
  def ReprocessPredictions(self, realigned_contigs): #{
    LogMsg(self, "Reprocessing %s predictions..." % self.description)
    self.num_predictions = 0
    self.write_pretty = True
    self.write_special = True
    #self.output_file.Close()
    #parser = CandidateGroupParserCls(self.output_file.path)
    #output_path = GetFilePath(self.options.output_dir, self.options.lib,
    #  self.ext)
    #self.output_file = FileBoxCls(output_path, "w",
    #  "cannot create %s output file" % self.description)
    self.outfiles["main"].Close()
    parser = CandidateGroupParserCls(self.outfiles["main"].path)
    self.AddOutFile("main", self.ext, self.description)
    for group in parser: #{
      DebugMsg(self, "Reprocessing group %i:" % group.id)
      event = CandidateGroupCls(group.DataString())
      for (index, member) in enumerate(group.members): #{
        DebugMsg(self, "%s %s %s(%ibp):%.2f %s(%ibp):%.2f %s" % (
          member.IDString(), member.contig_info.ToString(with_kform=False),
          member.align_info_A.ContigCoordsString(),
          member.align_info_A.ContigSpan(), member.align_info_A.identity,
          member.align_info_B.ContigCoordsString(),
          member.align_info_B.ContigSpan(), member.align_info_B.identity,
          member.MainGenesString()))
        realigned_contig = realigned_contigs[member.contig_info.id]
        #contig_seq = None
        #if (self.store_seq): #{
        #  contig_seq = contig_sequences[member.contig_info.id]
        #} end if
        #if (self.ReprocessMember(member, contig_alignment, contig_seq)): #{
        # Contigs that are intronic or intragenic will not have any
        #   contig-to-transcript alignments
        #if (None == realigned_contig.best_align): #{
        #  LogMsg(self, "  No alignments for contig %s" % member.contig_info.id)
        #  continue
        #} end if
        if (self.ReprocessMember(member, realigned_contig)): #{
          event.members.append(group.members[index])
        #} end if
      #} end for
      if (0 == len(event.members)): #{
        DebugMsg(self, "Group %i has no good members after realignment" %
          group.id)
        continue
      #} end if
      event.UpdateMemberCount()
      self.OutputEvent(event)
    #} end for
  #} end def
# end class
