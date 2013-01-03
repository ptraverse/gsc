#! /usr/bin/env python
"""
fastq.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
from utils.files_paths import FileBoxCls

# CONSTANTS

class FastqFileCls: #{
  def __init__(self, path, fail_msg="cannot open fastq file",
      log_info=None, line_delim="", maintain_case=False): #{
    self.file = FileBoxCls(path, "r", fail_msg)
    self.line_delim = line_delim
    self.curr_line = None
    self.log_info = log_info
    self.finished = False
    self.maintain_case = maintain_case
  #} end def

  def __del__(self): #{
    self.file.Close()
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  def next(self): #{
    if (self.finished): #{
      raise StopIteration
    #} end if
    new_seq = None
    try: #{
      if (None == self.curr_line): #{
        self.curr_line = self.file.next()
      #} end if
      if (not self.curr_line.startswith("@")): #{
        raise FastqError("improperly formatted fastq file: sequence id line "
          "must begin with \"@\": \"%s\"." % self.curr_line)
      #} end if
      if (" " in self.curr_line): #{
        (seq_id, seq_extra) = self.curr_line.lstrip("@").split(" ", 1)
      else:
        seq_id = self.curr_line.lstrip("@")
        seq_extra = None
      #} end if
      new_seq = SequenceAndQualityCls(seq_id, seq_extra)
      self.curr_line = self.file.next()
      while (not self.curr_line.startswith("+")): #{
        if ("" != new_seq.sequence): #{
          new_seq.sequence += self.line_delim
        #} end if
        new_seq.sequence += self.curr_line
        try: #{
          self.curr_line = self.file.next()
        except StopIteration:
          self.finished = True
          break
        #} end try
      #} end while
      #self.curr_line = self.file.next()
      if (not self.curr_line.startswith("+")): #{
        raise FastqError("improperly formatted fastq file: sequence must be "
          "separated from quality with \"+\", not \"%s\"." % self.curr_line)
      #} end if
      self.curr_line = self.file.next()
      while (not self.curr_line.startswith("@")): #{
        if ("" != new_seq.quality): #{
          new_seq.quality += self.line_delim
        #} end if
        new_seq.quality += self.curr_line
        try: #{
          self.curr_line = self.file.next()
        except StopIteration:
          self.finished = True
          break
        #} end try
      #} end while
      if (not self.maintain_case): #{
        new_seq.sequence = new_seq.sequence.upper()
      #} end if
      return new_seq
    except StopIteration, e:
      self.finished = True
      raise e
    #} end except
  #} end def

  def Reset(self): #{
    self.file.close()
    self.curr_line = None
    self.finished = False
  #} end def

  def Close(self): #{
    self.file.Close()
  #} end def
#} end class

class SequenceAndQualityCls: #{
  def __init__(self, id, extra): #{
    self.id = id
    self.extra = extra
    self.sequence = ""
    self.quality = ""
  #} end def

  def __len__(self): #{
    return len(self.sequence)
  #} end def

  def Output(self, maintain_case=False): #{
    output_str = "@%s" % self.id
    if (None != self.extra): #{
      output_str += " %s" % self.extra
    #} end if
    if (maintain_case): #{
      output_str += "\n%s" % self.sequence
    else:
      output_str += "\n%s" % self.sequence.upper()
    #} end if
    output_str += "\n+\n%s" % self.quality
    return output_str
  #} end def
#} end def

#### EXCEPTION CLASSES ####
class FastqError(MyError): #{
  pass
#} end class
