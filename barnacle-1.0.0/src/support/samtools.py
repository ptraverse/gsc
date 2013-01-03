#! /usr/bin/env python
"""
samtools.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os, time

# import custom modules
from utils.log import CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, CheckConfigCommands,
  GetCommand, NormalizeChrID)
from utils.messages import LogMsg, ExtremeDebugMsg
from utils.files_paths import FileBoxCls, CleanLine
from utils.subprocesses import (RunCommandFromString, RunPipeFromStrings,
  STRING_OUT, STREAM_OUT)
from SAM_record import SAM_record

class SAMToolsCls: #{
  def __init__(self, align_path, options, log_info=None,
      out_file=None, sort=False, paired=False): #{
    SetupMainClass(self, options, log_info=log_info)
    CheckConfigCommands(self, "samtools")
    if (hasattr(options, "output_dir")): #{
      self.output_dir = options.output_dir
      if (None == out_file): #{
        out_file = "sam_out_tmp"
      #} end if
      self.out_file_path = os.path.join(options.output_dir, out_file)
      self.err_file_path = "%s.err" % self.out_file_path
    #} end if
    self.align_path = align_path
    # filter unmapped query (0x0004) and qc-failed (0x0200) and
    #   PCR duplicates (0x0400) -- need to bitwise OR the flags
    #   0x0004 | 0x0200 | 0x0400 == 1540
    filter_flag = 0x0004 | 0x0200 | 0x0400
    if (paired): #{
      # also filter unmapped mate (0x0008)
      #   0x0004 | 0x0200 | 0x0400 | 0x0008 == 1548
      filter_flag = filter_flag | 0x0008
    else:
      # also filter secondary alignments (0x0100)
      #   0x0004 | 0x0200 | 0x0400 | 0x0100 == 1796
      filter_flag = filter_flag | 0x0100
    #} end if
    self.cmd_template = "%s view -F %i %s \"REGION\"" % (
        GetCommand(self, "samtools"), filter_flag, self.align_path)
    self.sort = sort
    self.use_chr = None
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def ShouldChromUseChr(self): #{
    samtools_cmd = "%s view %s -H" % (GetCommand(self, "samtools"),
      self.align_path)
    bam_header = RunCommandFromString(samtools_cmd, stdout=STRING_OUT,
      dpt=self.options.dpt)[0]
    header_lines = bam_header.split("\n")
    ExtremeDebugMsg(self, "Bam Header Lines:")
    for line in header_lines: #{
      ExtremeDebugMsg(self, CleanLine(line))
      if (line.startswith(r"@SQ")): #{
        if ("\t" in line): #{
          line_fields = line.split("\t")
        else:
          line_fields = line.split(" ")
        #} end if
        for line_field in line_fields: #{
          ExtremeDebugMsg(self, "Field: %s" % line_field)
          if (line_field.startswith("SN:")): #{
            return ("chr" in line_field)
          #} end if
        #} end for
      #} end if
    #} end for
    raise SAMToolsError("could not get any chromosome IDs from bam header")
  #} end def

  def DetermineChrUse(self): #{
    self.use_chr = self.ShouldChromUseChr()
  #} end def

  def GetCommand(self, region): #{
    if (None == self.use_chr):
      good_region = region
    else:
      good_region = NormalizeChrID(region, self.use_chr)
    #} end if
    return self.cmd_template.replace("REGION", good_region)
  #} end def

  def Run(self, region, retry_num=0): #{
    start_time = time.time()
    #command = self.cmd_template.replace("REGION", region)
    command = self.GetCommand(region)
    ExtremeDebugMsg(self, "Running samtools view:\n%s" % command)
    # open output and error files
    fail_msg = "cannot open temporary samtools output file"
    out_file = {'path': self.out_file_path, 'mode': "w", 'fail_msg': fail_msg}
    fail_msg = "cannot open temporary samtools error file"
    err_file = {'path': self.err_file_path, 'mode': "w", 'fail_msg': fail_msg}
    # run command
    retry = False
    try: #{
      if (self.sort): #{
        # set the LC_ALL environmental variable to ensure sort uses byte values
        os.environ["LC_ALL"] = "C"
        commands = [command, "sort --temporary-directory %s" % self.output_dir]
        #ExtremeDebugMsg(self, "PIPE: %s" % " | ".join(commands))
        retcode = RunPipeFromStrings(commands, stdout=out_file,
          stderr=err_file, dpt=self.options.dpt)
      else:
        retcode = RunCommandFromString(command, stdout=out_file,
          stderr=err_file, dpt=self.options.dpt)
      #} end if
    except OSError, e:
      LogMsg(self, "Error executing samtools: %s" % e)
      retry = True
      retcode = 0
    #} end try
    ExtremeDebugMsg(self, "Time spent running samtools: %s" %
      TimeSpent(start_time))
    if (0 > retcode): #{
      LogMsg(self, "Samtools was terminated by signal %i" % -retcode)
      retry = True
    elif (0 < retcode):
      LogMsg(self, "WARNING: samtools exited with exit code: %i" % retcode)
    #} end if
    if (retry): #{
      if (self.options.max_sam_retries > retry_num): #{
        LogMsg(self, "WARNING: Could not start samtools, will try again. " +
          "Retry Attempt: %i" % (retry_num + 1))
        # wait for ten minutes
        time.sleep(600);
        # try again
        self.Run(region, retry_num + 1)
        return
      else:
        self.PrintErrors()
        raise SAMToolsError \
          ("ERROR: could not run samtools command: %s" % command)
      #} end if
    #} end if
    if (self.PrintErrors()): #{
      raise SAMToolsError \
        ("ERROR: errors encountered while running samtools " +
         "command: %s" % command)
    #} end if
  #} end def

  def GetReads(self, region): #{
    #command = self.cmd_template.replace("REGION", region)
    command = self.GetCommand(region)
    ExtremeDebugMsg(self, "COMMAND: %s" % command)
    fail_msg = "cannot open temporary samtools error file"
    err_file = {'path': self.err_file_path, 'mode': "w", 'fail_msg': fail_msg}
    return RunCommandFromString(command, stdout=STREAM_OUT, stderr=err_file,
      dpt=self.options.dpt)
  #} end def

  def GetParsedReads(self, region, require_perfect=True,
      allow_softclipping=True, max_edit_distance=0): #{
    for raw_record in self.GetReads(region): #{
      raw_record = CleanLine(raw_record)
      #ExtremeDebugMsg(self, "%s" % raw_record)
      # use Rod's SAM_record code to parse reads from r2c alignment file
      record = SAM_record(raw_record)
      # make sure that the line was not a header (just in case)
      if (record.header): #{
        ExtremeDebugMsg(self, "Line \"%s\" is a header line, it will be "
          "skipped..." % record.raw_sam)
        continue
      #} end if
      # make sure the read is mapped (just in case)
      if (record.q_unmapped()): #{
        ExtremeDebugMsg(self, "Read %s is unmapped, will be skipped..." %
          record.qname)
        continue
      #} end if
      # make sure the read is a perfect match, possibly with soft-clipping
      if (require_perfect): #{
        ExtremeDebugMsg(self, "Requiring perfect match! Allow soft-clipping: "
            "%s" % allow_softclipping)
        if (allow_softclipping): #{
          if (not record.perfect_softclipped(max_edit_distance)): #{
            ExtremeDebugMsg(self, "Skipping imperfect soft-clipped alignment "
              "for read %s: %s ED:%i" % (record.qname, record.cigar,
                record.edit_distance))
            continue
          #} end if
        elif (not record.perfect_match(max_edit_distance)): #{
          ExtremeDebugMsg(self, "Skipping imperfect alignment for read %s: %s "
            "ED:%i" % (record.qname, record.cigar, record.edit_distance))
          continue
        #} end if
      #} end if
      #if (None != max_edit_distance and
      #    max_edit_distance < record.edit_distance): #{
      #  #ExtremeDebugMsg(self, "Skipping alignment for read %s with edit "
      #  #  "distance %i (max %i)" % (record.qname, record.edit_distance,
      #  #  max_edit_distance))
      #  continue
      #} end if
      #ExtremeDebugMsg(self, "%s" % record)
      yield record
    #} end for
  #} end def

  def PrintErrors(self): #{
    errors = False
    fail_msg = "could not open samtools error file"
    err_file = FileBoxCls(self.err_file_path, "r", fail_msg)
    for line in err_file: #{
      LogMsg(self, line)
      errors = True
    #} end for
    err_file.close()
    return errors
  #} end def
#} end class

def DecodeFlag(flag_value): #{
  flag_fields = {0x0001: "read is paired in sequencing",
    0x0002: "read is mapped in a proper pair", 0x0004: "unmapped query",
    0x0008: "unmapped mate", 0x0010: "query on reverse strand",
    0x0020: "mate on reverse strand", 0x0040: "first read",
    0x0080: "second read", 0x0100: "not primary",
    0x0200: "fails qc", 0x0400: "duplicate read",
  }
  for field in sorted(flag_fields.keys()): #{
    if (field & flag_value): #{
      print flag_fields[field]
    #} end if
  #} end for
#} end def

class SAMToolsError(MyError): #{
  pass
#} end class
