#! /usr/bin/env python
"""
get_primer3_input.py

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
from utils.general import SetupMainClass, TimeSpent, WriteCommand
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls)
from parsers.candidate_group_parser import CandidateGroupParserCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

class PrimerDesignerCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    # self.members[ctg_id] = member_list
    self.members = dict()
    # self.contigs[ctg_id] = sequence
    self.contigs = dict()
    # self.transcripts[gene_id] = sequence
    self.transcripts = dict()
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def DesignPrimers(self): #{
    LogMsg(self, "Desiging primers...")
    start = time.time()
    # load the data needed to design primers
    self.LoadData()
    # create the output file
    output_path = os.path.join(self.options.output_dir, "%s.primers" %
      os.path.basename(self.options.barnacle_path))
    DebugMsg("Writing results to: %s" % output_path)
    self.primers_file = FileBoxCls(output_path, "w",
      "cannot create primer output file")
    for member_list in self.members.itervalues(): #{
      for member in member_list: #{
        self.DesignPrimerForMember(member)
      #} end for
    #} end for
    LogMsg(self, "Time spent desiging primers: %s" % TimeSpent(start))
  #} end def

  def LoadData(self): #{
    LogMsg(self, "Loading data...")
    start = time.time()
    # load the prediction data
    LogMsg(self, "Loading Barnacle predictions...")
    group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    for group in group_parser: #{
      if ("fusion" != group.event_type): #{
        raise PrimerDesignerError("unsupported event type: %s" %
          group.event_type)
      #} end if
      member = self.SelectMember(group)
      DebugMsg(self, "Selected member:\n%s" % member.DataString())
      self.PreProcessMember(member, group.event_type)
    #} end for
    # load the contig sequences
    self.LoadContigSequences()
    # load the alignments for the chosen members
    self.LoadAlignments()
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
#} end class

class PrimerMemberCls: #{
  pass
#} end class

#### EXCEPTION CLASSES ####
class PrimerDesignerError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Use Primer3 and gfPcr to create qPCR and sequence "
    "validation primers for the Barnacle events in the given file.")
  args = [ "LIB", "BARNACLE_FILE", "ALIGNMENTS_FILES", "CONTIG_SEQS_FILE",
    "TRANSCRIPT_SEQS_FILE", "TRANSCRIPTOME_HOSTNAME", "TRANSCRIPTOME_PORT",
    "GENOME_HOSTNAME", "GENOME_PORT" ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--qpcr-buffer",
                    type="int", metavar="N",
                    help="Ensure that qPCR primers are at least Nbp from the "
                      "breakpoint. [default: %default]")
  parser.add_option("--seq-buffer",
                    type="int", metavar="N",
                    help="Ensure that sequencing primers are at least Nbp "
                      "from the breakpoint. [default: %default]")
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(qpcr_buffer=20,
                      seq_buffer=70,
                      force=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  try: #{
    int(options.t_port)
  except ValueError, e:
    ErrMsg("  Invalid transcriptome port: \"%s\", " % options.t_port +
      "must be a positive integer.")
    opts_good = False
  #} end if
  try: #{
    int(options.g_port)
  except ValueError, e:
    ErrMsg("  Invalid genome port: \"%s\", " % options.g_port +
      "must be a positive integer.")
    opts_good = False
  #} end if
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
    "primers")
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "primers", options.output_dir)
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
    options.t_host = args[5]
    options.t_port = args[6]
    options.g_host = args[7]
    options.g_port = args[8]
    if (CheckPaths(options)): #{
      try: #{
        primer_designer = PrimerDesignerCls(options)
        WriteCommand(primer_designer, sys.argv)
        primer_designer.DesignPrimers()
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
      "alignment files (ALIGNMENTS_FILES); the path to a fasta-formatted "
      "candidate contig sequence file (CONTIG_SEQS_FILE); the path to a "
      "fasta-formatted transcript sequence file (TRANSCRIPT_SEQS_FILE); the "
      "hostname and port of a server for aligning potential primers to the "
      "transcriptome (TRANSCRIPTOME_HOSTNAME, TRANSCRIPTOME_PORT); and the "
      "hostname and port of a server for aligning potential primers to the "
      "genome (GENOME_HOSTNAME, GENOME_PORT)")
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
