#! /usr/bin/env python
"""
candidate_group_parser.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import os

# import custom modules
from utils.error import MyError
from utils.files_paths import FileBoxCls, GetFilePath, CheckFilePath
from common.candidate_group import GroupParserCls

# HOW TO USE THIS CLASS
# create an object:
#   parser = CandidateGroupParserCls("path_to_data_file")
# 1) load all groups into memory
#   groups = parser.ParseDataFile()
#   for group in groups:
#     foo(group)
#   # end for
# OR
# 2) load only one group into memory at a time
#   for group in parser:
#     foo(group)
#   # end for

# CONSTANTS
TOPOLOGY_EXT_LIST = list([
  ("unknown",                            "unk"),
  ("end-duplication",                    "edu"),
  ("intrachr-non-colinear",              "ncl"),
  ("junction-duplication",               "jdu"),
  ("intrachr-same-strand",               "sam"),
  ("read-through",                       "rth"),
  ("intrachr-opp-strand",                "opp"),
  ("local-inversion",                    "lin"),
  ("interchr",                           "dif"),
  ("multi",                              "mul"),
  ("gap-tandem-duplication",             "gdu"),
  ("gap-nontandem-duplication",          "gnd"),
  ("gap-tandem-inverted_duplication",    "gid"),
  ("gap-nontandem-inverted_duplication", "gnid"),
  ("gap-intronic_BLAT_miss",             "gbm"),
  ("gap-internal_inversion",             "gii"),
  ("gap-genic_rearrangement",            "ggr"),
  ("gap-genic_inversion",                "ggi"),
])

class CandidateGroupParserCls: #{
  def __init__(self, data_file_path, keep_lines=False, check_data=False): #{
    CheckFilePath(data_file_path, "candidate group file")
    self.group_parser = GroupParserCls(keep_lines=keep_lines)
    self.check_data = check_data
    fail_message = "cannot open data file"
    self.data_file = FileBoxCls(data_file_path, "r", fail_message)
    self.groups = list()
  #} end def

  def __del__(self): #{
    # close data file if it is open
    self.CloseDataFile()
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  # Load the entire data file into memory
  # Do not mix with using GetNextGroup() method
  def ParseDataFile(self): #{
    #self.OpenDataFile()
    for group_line in self.data_file: #{
      #group_line = CleanLine(group_line)
      # skip blank lines
      #if ("" == group_line): #{
      #  continue
      #} end if
      self.groups.append(self.group_parser.ParseGroup \
        (group_line, self.data_file, check_data=self.check_data))
    #} end for
    self.CloseDataFile()
    return self.groups
  #} end def

  # Load a single group from the data file into memory
  # Do not mix with using ParseDataFile() method
  def GetNextGroup(self): #{
    return self.next()
  #} end def

  def next(self): #{
    #if (None == self.data_file): #{
    #  self.OpenDataFile()
    #} end if
    group_line = ""
    # skip blank lines
    while ("" == group_line): #{
      #group_line = CleanLine(self.data_file.next())
      group_line = self.data_file.next()
    #} end if
    return self.group_parser.ParseGroup \
      (group_line, self.data_file, check_data=self.check_data)
  #} end def

  def Close(self): #{
    self.CloseDataFile()
  #} end def

  def CloseDataFile(self): #{
    if (not hasattr(self, "data_file")): #{
      return
    #} end if
    if (None == self.data_file): #{
      return
    #} end if
    if (self.data_file.closed): #{
      return
    #} end if
    self.data_file.close()
    #self.data_file = None
  #} end def

  def close(self): #{
    self.CloseDataFile()
  #} end def

  def GroupLine(self): #{
    if (not self.group_parser.keep_lines): #{
      raise CandidateGroupParserError \
        ("cannot get group line when keep_lines flag was not set")
    #} end if
    return self.group_parser.group_line
  #} end def

  def MemberLines(self): #{
    if (not self.group_parser.keep_lines): #{
      raise CandidateGroupParserError \
        ("cannot get member lines when keep_lines flag was not set")
    #} end if
    return self.group_parser.member_lines
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class CandidateGroupParserError(MyError): #{
  pass
#} end class
