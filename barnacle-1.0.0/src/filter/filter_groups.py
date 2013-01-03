#! /usr/bin/env python
"""
filter_groups.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import math, os, re, sys, time, traceback

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
from utils.general import (SetupMainClass, TimeSpent, WriteCommand)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls)
from common.candidate_group import CandidateGroupCls
from parsers.candidate_group_parser import (TOPOLOGY_EXT_LIST,
                                            CandidateGroupParserCls)

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "FILTER SUCCESS"
MSG_FAIL = "FILTER FAIL"

# Filter types
FT_BOOL    = "bool"
FT_INT     = "int"
FT_PERCENT = "percent"
FT_FLOAT   = "float"
FT_SET     = "set"

# Filter descriptions
FD_READ_TO_CTG    = ("Minimum number of supporting read to contig alignments")
FD_READ_TO_CTG_U  = ("Minimum number of unique supporting read to contig "
                     "alignments")
FD_READ_TO_CTG_RU = ("Require unique supporting read to contig alignments")
FD_PAIR_TO_GEN    = ("Minimum number of supporting read-pair to genome "
                     "alignments")
FD_PAIR_TO_GEN_F  = ("Minimum number of mapq-filtered supporting read-pair "
                     "to genome alignments")
FD_PAIR_TO_GEN_I  = ("Discard intronic supporting read-pair to genome "
                     "alignments")
FD_RNA            = ("Fail RNA")
FD_MITO           = ("Fail mitochondrial DNA")
FD_NUM_GROUPS     = ("Maximum number of groups per contig")
FD_ALIGN          = ("Require multiple aligners")
FD_PID            = ("Minimum percent identity")
FD_ADJUST_PID     = ("Adjust minimum percent identity to be less stringent "
                     "for shorter alignments")
FD_MISMATCHES     = ("Number of mismatches to use when adjusting minimum "
                     "percent identity")
FD_RUN            = ("Fail single base run gap events")
FD_RUN_SOFT       = ("Allow one divergent base in runs")
FD_CTG_OLAP       = ("Maximum overlap in contig of alignments for "
                     "paired-alignment candidates.")
FD_CTG_REP        = ("Minimum fraction of contig length represented by "
                     "alignments")
FD_MULTI_MAPS     = ("Fail candidates with no unique contig-to-genome "
                     "alignments")
FD_REPEATS        = ("Fail groups overlapping annotated repeat regions")
FD_POLY_A         = ("Fail groups that are probably polyA tails "
                     "(single-base runs at the very beginning or end of a "
                     "contig)")
FD_IG_TOPOLOGIES  = ("Alignment topologies to fail")
FD_NO_GOOD_MEMBS  = ("Fail groups with no members passing all filters")

# sort methods
SORT_BY_R2CU = 0

class GroupFilteringCls: #{
  def __init__(self, options, log_info=None): #{
    SetupMainClass(self, options, log_info=log_info)
    self.num_groups = 0
    self.num_passing = 0
    self.num_warnings = 0
    self.output_files = dict()
    self.SelectSortMethod(options)
    self.CreateFilterDictionary(options)
    self.group_parser = CandidateGroupParserCls(self.options.barnacle_path)
    self.CompileRegexes()
    #self.gene_names = None
  #} end def

  def __del__(self): #{
    for file_box in self.output_files.itervalues(): #{
      file_box.Close()
    #} end for
    CloseLogFile(self)
  #} end def

  def SelectSortMethod(self, options): #{
    self.sort_method = None
    if (options.sort_by_r2cu): #{
      self.sort_method = SORT_BY_R2CU
      return
    #} end if
  #} end def

  def CreateFilterDictionary(self, options): #{
    LogMsg(self, "Creating filters...")
    self.filters = dict()
    self.AddFilter("READ_TO_CTG", FD_READ_TO_CTG, "READS_TO_CONTIG",
                   FT_INT, "read_to_ctg")
    self.AddFilter("READ_TO_CTG_U", FD_READ_TO_CTG_U, "READS_TO_CONTIG",
                   FT_INT, "read_to_ctg_u")
    self.AddFilter("READ_TO_CTG_U_REQ", FD_READ_TO_CTG_RU, "READS_TO_CONTIG",
                   FT_BOOL, "req_u_read_to_ctg")
    self.AddFilter("PAIR_TO_GEN", FD_PAIR_TO_GEN, "PAIR_SUPPORT",
                   FT_INT, "pair_to_gen")
    self.AddFilter("PAIR_TO_GEN_F", FD_PAIR_TO_GEN_F, "PAIR_SUPPORT",
                   FT_INT, "pair_to_gen_f")
    self.AddFilter("PAIR_TO_GEN_I", FD_PAIR_TO_GEN_I, "PAIR_SUPPORT",
                   FT_BOOL, "ignore_intron")
    self.AddFilter("RNA", FD_RNA, "RNA_REPEAT",
                   FT_BOOL, "no_rna")
    self.AddFilter("MITO", FD_MITO, "MITOCHONDRIAL",
                   FT_BOOL, "no_mito")
    self.AddFilter("NUM_GROUPS", FD_NUM_GROUPS, "TOO_MANY_GROUPS",
                   FT_INT, "max_num_groups")
    self.AddFilter("ALIGN", FD_ALIGN, "SINGLE_ALIGN",
                   FT_BOOL, "reguire_multiple_aligns")
    self.AddFilter("PID", FD_PID, "PERCENT_IDENTITY",
                   FT_PERCENT, "min_pid")
    self.AddFilter("ADJUST_PID", FD_ADJUST_PID, "PERCENT_IDENTITY",
                   FT_BOOL, "adjust_pid")
    self.AddFilter("MISMATCHES", FD_MISMATCHES, "PERCENT_IDENTITY",
                   FT_INT, "max_mismatches")
    self.AddFilter("RUN", FD_RUN, "GAP_RUN",
                   FT_BOOL, "no_runs")
    self.AddFilter("RUN_SOFT", FD_RUN_SOFT, "GAP_RUN",
                   FT_BOOL, "soft_runs")
    self.AddFilter("CTG_OLAP", FD_CTG_OLAP, "CONTIG_OVERLAP",
                   FT_INT, "ctg_olap")
    self.AddFilter("CTG_REP", FD_CTG_REP, "CONTIG_REPRESENTED",
                   FT_FLOAT, "ctg_rep")
    self.AddFilter("MULTI_MAPS", FD_MULTI_MAPS, "MULTI_MAPS",
                   FT_BOOL, "filter_multi_maps")
    self.AddFilter("REPEATS", FD_REPEATS, "REPEAT_SEQ",
                   FT_BOOL, "filter_repeats")
    self.AddFilter("POLY_A", FD_POLY_A, "POLY_A",
                   FT_BOOL, "filter_polyA")
    self.AddFilter("IG_TOPOLOGIES", FD_IG_TOPOLOGIES, "IGNORED_TOPOLOGY",
                   FT_SET, "ignore_topologies")
    self.AddFilter("NO_GOOD_MEMBERS", FD_NO_GOOD_MEMBS, "NO_GOOD_MEMBERS",
                   FT_BOOL, "member_wise")
  #} end def

  def AddFilter(self, key, description, fail_reason, filter_type,
      attr_name): #{
    try: #{
      value = getattr(self.options, attr_name)
    except AttributeError:
      value = None
    #} end try
    if ("IG_TOPOLOGIES" == key and None != value): #{
      if ("" == value or "none" == value.lower()): #{
        value = None
      else:
        value = set(value.split(","))
      #} end if
    #} end if
    self.filters[key] = FilterCls(description, fail_reason, filter_type, value)
  #} end def

  def CompileRegexes(self): #{
    if (self.filters['RUN'].value): #{
      self.run_regex_hard = re.compile(r"^(.)\1*$")
      self.run_regex_soft1 = re.compile(r"^(.)\1*.\1*$")
      self.run_regex_soft2 = re.compile(r"^.(.)\1*$")
    #} end if
  #} end def

  def FilterFile(self, print_filters=True): #{
    start_time = time.time()
    LogMsg(self, "Filtering %s..." % self.options.barnacle_path)
    # get the "acceptable" gene names
    #self.GetRefGeneNames()
    self.SetupOutputFiles()
    if (print_filters): #{
      self.PrintFilters()
    #} end if
    if (None != self.sort_method): #{
      groups = list()
    #} end if
    for group in self.group_parser: #{
      #DebugMsg(self, group.FullDataString(), newline=False)
      self.FilterGroup(group)
      #group.FilterGeneNames(self.gene_names)
      if (None == self.sort_method): #{
        self.OutputGroup(group)
      else:
        groups.append(group)
      #} end if
      if (self.options.check_data): #{
        if (0 < len(group.member_warnings)): #{
          for warning in group.member_warnings: #{
            LogMsg(self, "WARNING: %s" % warning)
          #} end for
          self.num_warnings += len(group.member_warnings)
        #} end if
      #} end if
      self.num_groups += 1
      if (group.PassedFilters()): #{
        self.num_passing += 1
      #} end if
      DebugMsg(self, "*"*80)
    #} end for
    self.group_parser.Close()
    LogMsg(self, "Time spent filtering groups: %s" % TimeSpent(start_time))
    if (None != self.sort_method): #{
      self.SortGroups(groups)
      self.OutputGroups(groups)
    #} end if
  #} end def

  #def GetRefGeneNames(self): #{
  #  self.gene_names = GetGeneNamesFromFile(self.options.gene_names_path,
  #    self.log_info)
  #} end def

  def SetupOutputFiles(self): #{
    output_file_list = [("passing", "pass")]
    if (self.options.print_fails): #{
      output_file_list.append(("failing", "fail"))
    #} end if
    if (self.options.split_out): #{
      output_file_list.extend(TOPOLOGY_EXT_LIST)
    #} end if
    for output_type, ext in output_file_list: #{
      self.output_files[output_type] = FileBoxCls(self.OutputPath(ext), "w",
        "could not open output file for %s groups" % output_type)
    #} end for
  #} end def

  def PrintFilters(self): #{
    filter_file = FileBoxCls(self.OutputPath("filters"), "w",
      "could not open filter file")
    try: #{
      for filter_name in sorted(self.filters.keys()): #{
        #LogMsg(self, "%s: %s" %
        #  (self.filters[filter_name].description,
        #   self.filters[filter_name].ValueString()))
        filter_file.WriteLine("%s: %s" %
          (self.filters[filter_name].description,
           self.filters[filter_name].ValueString()))
      #} end for
    finally:
      filter_file.close()
    #} end try
  #} end def

  def DisplayFilters(self): #{
    for filter_name in sorted(self.filters.keys()): #{
      LogMsg(self, "%s: %s" %
        (self.filters[filter_name].description,
         self.filters[filter_name].ValueString()))
    #} end for
  #} end def

  def FilterGroup(self, group): #{
    DebugMsg(self, "*"*80 + "\nFiltering group: %i..." % group.id)
    # clear the fail reasons for the group
    group.fail_reasons = set()
    # set filters that pass as long as one member passes to fail by default
    self.SetPassByOneFailures(group)
    # apply group-level filters
    self.ApplyAlignerFilter(group)
    self.ApplyMitoFilter(group)
    self.ApplyPairToGenomeFilter(group)
    self.ApplyRNAFilter(group)
    self.ApplyRepeatFilter(group)
    self.ApplyTopologyFilter(group)
    # check filters that pass as long as one member passes, removing fail
    # reasons as passing member(s) are found
    num_good_members = 0
    DebugMsg(self, "-"*80)
    for member in group.members: #{
      if (group.PassedFilters() and not self.options.member_wise): #{
        break
      #} end if
      self.ApplyContigOverlapFilter(group, member)
      self.ApplyNumGroupsFilter(group, member)
      self.ApplyIdentityFilter(group, member)
      self.ApplyReadToContigFilter(group, member)
      self.ApplyRunFilter(group, member)
      self.ApplyContigRepresentedFilter(group, member)
      self.ApplyMultiMapFilter(group, member)
      self.ApplyPolyAFilter(group, member)
      self.ApplyMemberTopologyFilter(group, member)
      if (self.options.member_wise and member.PassedFilters()): #{
        num_good_members += 1
      #} end if
      DebugMsg(self, "Member %s: %s\n" %
        (member.IDString(), member.Status()) + "-"*80)
    #} end for
    if (self.options.member_wise): #{
      if (0 == num_good_members): #{
        self.AddFailReason(group, 'NO_GOOD_MEMBERS')
      #} end if
      #if (None != self.filters['IG_TOPOLOGIES'].value and
      #    group.PassedFilters() and
      #    num_good_members < group.num_members):
      #  group.UpdateTopologies()
      #} end if
      #if (group.PassedFilters()): #{
      #  group.UpdateMemberCount(num_good_members)
      #} end if
    #} end if
    DebugMsg(self, "Group Status: %s" % group.Status())
  #} end def

  def SetPassByOneFailures(self, group): #{
    if (0 < self.filters['NUM_GROUPS'].value): #{
      self.AddFailReason(group, 'NUM_GROUPS')
    #} end if
    self.AddFailReason(group, 'PID')
    self.AddFailReason(group, 'READ_TO_CTG')
    self.AddFailReason(group, 'CTG_REP')
    if (self.filters['MULTI_MAPS'].value): #{
      self.AddFailReason(group, 'MULTI_MAPS')
    #} end if
    if (group.gap): #{
      if (self.filters['RUN'].value): #{
        self.AddFailReason(group, 'RUN')
      #} end if
      if (self.filters['POLY_A'].value): #{
        self.AddFailReason(group, 'POLY_A')
      #} end if
    else:
      self.AddFailReason(group, 'CTG_OLAP')
    #} end if
  #} end def

  def ApplyAlignerFilter(self, group): #{
    # gapped events are only looked for in blat alignments
    if (self.filters['ALIGN'].value and
        2 > group.num_aligners      and
        not group.gap):
      self.AddFailReason(group, 'ALIGN')
    #} end if
  #} end def

  def ApplyMitoFilter(self, group): #{
    if (self.filters['MITO'].value): #{
      chromosomes = group.coord_pair_A.chrom
      if (None != group.coord_pair_B): #{
        chromosomes += group.coord_pair_B.chrom
      #} end if
      if ("M" in chromosomes): #{
        self.AddFailReason(group, 'MITO')
      #} end if
    #} end if
  #} end def

  def ApplyPairToGenomeFilter(self, group): #{
    # if ignoring "intronic" use only "exonic" read-pair support
    if (self.filters['PAIR_TO_GEN_I'].value): #{
      p2g_all      = group.pair_to_genome_exonic
      p2g_filtered = group.pair_to_genome_filtered_exonic
    # otherwise use total read-pair support
    else:
      p2g_all      = group.pair_to_genome_all
      p2g_filtered = group.pair_to_genome_filtered_all
    #} end if
    if (self.filters['PAIR_TO_GEN'].value   > p2g_all or
        self.filters['PAIR_TO_GEN_F'].value > p2g_filtered):
      self.AddFailReason(group, 'PAIR_TO_GEN')
    #} end if
  #} end def

  def ApplyRNAFilter(self, group): #{
    if (self.filters['RNA'].value and group.rna): #{
      self.AddFailReason(group, 'RNA')
    #} end if
  #} end def

  def ApplyRepeatFilter(self, group): #{
    if (self.filters['REPEATS'].value): #{
      if (group.OverlapsRepeat()): #{
        self.AddFailReason(group, 'REPEATS')
      #} end if
    #} end if
  #} end def

  def ApplyTopologyFilter(self, group): #{
    if (None == self.filters['IG_TOPOLOGIES'].value): #{
      return
    #} end if
    group_topologies = set(group.topologies)
    ignored_topologies = self.filters['IG_TOPOLOGIES'].value
    if (group_topologies.issubset(ignored_topologies)): #{
      self.AddFailReason(group, 'IG_TOPOLOGIES')
    #} end if
  #} end def

  def ApplyMemberTopologyFilter(self, group, member): #{
    if (None == self.filters['IG_TOPOLOGIES'].value): #{
      return
    #} end if
    if (member.topology in self.filters['IG_TOPOLOGIES'].value): #{
      self.AddFailReason(member, 'IG_TOPOLOGIES')
    #} end if
  #} end if

  def ApplyNumGroupsFilter(self, group, member): #{
    if (self.filters['NUM_GROUPS'].value >= member.num_groups): #{
      self.RemoveFailReason(group, 'NUM_GROUPS')
    else:
      self.AddFailReason(member, 'NUM_GROUPS')
    #} end if
  #} end def

  def ApplyIdentityFilter(self, group, member): #{
    if (self.filters['ADJUST_PID'].value): #{
      min_pid_A = self.AdjustMinPID(member.align_info_A.ContigSpan())
      #if (member.gap and member.GapIsInternal()): #{}
      if (member.gap): #{
        min_pid_A = self.GapAdjustMinPID(min_pid_A, member)
      #} end if
      min_pid_B = self.AdjustMinPID(member.align_info_B.ContigSpan())
    else:
      min_pid_A = self.filters['PID'].value
      min_pid_B = self.filters['PID'].value
    #} end if
    DebugMsg(self, "Minimum PID: %.2f; %.2f" % (min_pid_A, min_pid_B))
    if (min_pid_A <= member.align_info_A.identity and
        min_pid_B <= member.align_info_B.identity):
      self.RemoveFailReason(group, 'PID')
    else:
      self.AddFailReason(member, 'PID')
    #} end if
  #} end if

  def AdjustMinPID(self, align_len): #{
    min_match = max(1, align_len-self.filters['MISMATCHES'].value)
    raw_adjusted_pid = (float(min_match) / float(align_len))
    adjusted_pid = math.floor(raw_adjusted_pid*1000)/10.0
    min_pid = min(adjusted_pid, self.filters['PID'].value)
    DebugMsg(self, "raw: %f, adjusted: %f, given: %f, min: %f" %
      (raw_adjusted_pid, adjusted_pid, self.filters['PID'].value, min_pid))
    return min_pid
  #} end def

  def GapAdjustMinPID(self, min_pid, member): #{
    qAliSize = member.align_info_A.ContigSpan()
    tAliSize = member.align_info_A.GenomeSpan()
    # only count the size difference if qAliSize is bigger
    sizeDif = max(0, qAliSize - tAliSize)
    # assume that the gap is the only internal portion of
    # the contig that is not aligned
    #total = qAliSize - member.GapSpan()
    # use the target blocks to get the number of bases aligned
    total = 0
    for block in member.blocks_A: #{
      (block_start, block_end) = block
      total += abs(block_start-block_end)+1
    #} end for
    # if the gap is at the edge, assume that there are no insertions
    insertFactor = 0
    # if the gap is internal, assume it is the only insertion
    if (member.GapIsInternal()): #{
      insertFactor = 1
    #} end if
    milliBad = float(1000 * (self.filters['MISMATCHES'].value + insertFactor +
      round(3*math.log(1+sizeDif)))) / float(total)
    raw_adjusted_pid = 100.0-(milliBad*0.1)
    adjusted_pid = math.floor(raw_adjusted_pid*10)/10.0
    min_pid = min(adjusted_pid, min_pid)
    DebugMsg(self, "gap-raw: %f, gap-adjusted: %f, given: %f, min: %f" %
      (raw_adjusted_pid, adjusted_pid, self.filters['PID'].value, min_pid))
    return min_pid
  #} end def

  def ApplyReadToContigFilter(self, group, member): #{
    if (self.filters['READ_TO_CTG_U'].value <= member.read_to_ctg_unique or
       (self.filters['READ_TO_CTG'].value   <= member.read_to_contig_all and
        not self.filters['READ_TO_CTG_U_REQ'].value)):
      self.RemoveFailReason(group, 'READ_TO_CTG')
    else:
      self.AddFailReason(member, 'READ_TO_CTG')
    #} end if
  #} end def

  def ApplyRunFilter(self, group, member): #{
    # only apply the run filter to gapped members
    if (not member.gap or not self.filters['RUN'].value): #{
      return
    #} end if
    sequence = member.event_seq.lower()
    if ((not self.filters['RUN_SOFT'].value and
         None == self.run_regex_hard.search(sequence)) or
        (None == self.run_regex_soft1.search(sequence) and
         None == self.run_regex_soft2.search(sequence))):
      #LogMsg(self, "Removing RUN filter failure. Seq %s" %
      #    member.event_seq)
      self.RemoveFailReason(group, 'RUN')
    else:
      self.AddFailReason(member, 'RUN')
    #} end if
  #} end def

  def ApplyContigOverlapFilter(self, group, member): #{
    if (member.gap): #{
      return
    #} end if
    if (self.filters['CTG_OLAP'].value >=
        member.meta_fields['ctg_overlap']): #{
      self.RemoveFailReason(group, 'CTG_OLAP')
    else:
      self.AddFailReason(member, 'CTG_OLAP')
    #} end if
  #} end def

  def ApplyContigRepresentedFilter(self, group, member): #{
    if (self.options.recalculate_ctg_rep): #{
      ctg_represented_fraction = AdHocGetContigRepresented(member,
        self.log_info, self.options.robust_crc)
    else:
      ctg_represented_fraction = member.TotalContigRepresented()
    #} end if
    if (self.filters['CTG_REP'].value <= ctg_represented_fraction): #{
      self.RemoveFailReason(group, 'CTG_REP')
    else:
      self.AddFailReason(member, 'CTG_REP')
    #} end def
  #} end def

  def ApplyMultiMapFilter(self, group, member): #{
    if (not self.filters['MULTI_MAPS'].value): #{
      return
    #} end if
    if (0 == member.MultiMapped()): #{
      self.RemoveFailReason(group, 'MULTI_MAPS')
    else:
      self.AddFailReason(member, 'MULTI_MAPS')
    #} end if
  #} end def

  def ApplyPolyAFilter(self, group, member): #{
    if (not member.gap or not self.filters['POLY_A'].value): #{
      return
    #} end if
    if (not member.PolyASuspect()): #{
      self.RemoveFailReason(group, 'POLY_A')
    else:
      self.AddFailReason(member, 'POLY_A')
    #} end if
  #} end def

  def AddFailReason(self, group, fail_key): #{
    group.fail_reasons.add(self.filters[fail_key].fail_reason)
  #} end def

  def RemoveFailReason(self, group, fail_key): #{
    group.fail_reasons.discard(self.filters[fail_key].fail_reason)
  #} end def

  def OutputGroup(self, src_group): #{
    pass_group = CandidateGroupCls(src_group.DataString())
    fail_group = CandidateGroupCls(src_group.DataString())
    for member in src_group.members: #{
      # add the member to the appropriate pass or fail group
      if (member.PassedFilters() and src_group.PassedFilters()): #{
        pass_group.members.append(member)
      else:
        fail_group.members.append(member)
      #} end if
    #} end for
    if (src_group.PassedFilters()): #{
      pass_group.UpdateMemberCount()
      if (None != self.filters['IG_TOPOLOGIES'].value and
          pass_group.num_members < pass_group.total_members): #{
        pass_group.UpdateTopologies()
      #} end if
      group_string = pass_group.FullDataString(readable=self.options.pretty)
      DebugMsg(self, "PASS Output:\n%s" % group_string, newline=False)
      self.output_files['passing'].Write(group_string)
      if (self.options.split_out): #{
        if (1 == len(group.topologies)): #{
          output_type = group.topologies[0]
        else:
          output_type = 'multi'
        #} end if
        self.output_files[output_type].Write(group_string)
      #} end if
    #} end if
    if (self.options.print_fails and 0 < len(fail_group.members)): #{
      fail_group.UpdateMemberCount()
      if (None != self.filters['IG_TOPOLOGIES'].value and
          fail_group.num_members < fail_group.total_members): #{
        fail_group.UpdateTopologies()
      #} end if
      group_string = fail_group.FullDataString(readable=self.options.pretty)
      DebugMsg(self, "FAIL Output:\n%s" % group_string, newline=False)
      self.output_files['failing'].Write(group_string)
    #} end if
  #} end def

  #def OutputGroup(self, group): #{
  #  if (self.options.pretty): #{
  #    group_data_list = [group.ReadableString()]
  #  else:
  #    group_data_list = [group.DataString()]
  #  #} end if
  #  for member in group.members: #{
  #    # only add the member if it passed or the group failed
  #    if (member.PassedFilters() or not group.PassedFilters()): #{
  #      if (self.options.pretty): #{
  #        group_data_list.append(member.ReadableString())
  #      else:
  #        group_data_list.append(member.DataString())
  #      #} end if
  #    #} end if
  #  #} end for
  #  group_data_string = "\n".join(group_data_list)+"\n\n"
  #  DebugMsg(self, "Output:\n%s" % group_data_string, newline=False)
  #  if (group.PassedFilters()): #{
  #    self.output_files['passing'].Write(group_data_string)
  #    if (self.options.split_out): #{
  #      if (1 == len(group.topologies)): #{
  #        output_type = group.topologies[0]
  #      else:
  #        output_type = 'multi'
  #      #} end if
  #      self.output_files[output_type].Write(group_data_string)
  #    #} end if
  #  elif (self.options.print_fails):
  #    self.output_files['failing'].Write(group_data_string)
  #  #} end if
  #} end def

  def DisplayFilterResults(self): #{
    LogMsg(self, "%i groups examined\n" % self.num_groups +
      "%i groups passed all filters" % self.num_passing)
    if (self.options.check_data): #{
      LogMsg(self, "%i warnings encountered" % self.num_warnings)
    #} end if
  #} end def

  def SortGroups(self, groups): #{
    if (SORT_BY_R2CU == self.sort_method): #{
      groups.sort(key=lambda group: group.max_read_to_ctg_u, reverse=True)
    else:
      raise GroupFilteringError("unrecognized sorting method: %i" %
        self.sort_method)
    #} end if
  #} end def

  def OutputGroups(self, groups): #{
    for group in groups: #{
      self.OutputGroup(group)
    #} end for
  #} end def

  def OutputPath(self, output_ext): #{
    input_file_name  = os.path.basename(self.options.barnacle_path)
    extension_recognized = False
    recognized_extensions = ["data", "pass", "fus", "fed", "ptd", "ped", "itd"]
    for input_ext in recognized_extensions: #{
      if (input_file_name.endswith(input_ext)): #{
        output_file_name = input_file_name.replace(input_ext, output_ext)
        extension_recognized = True
        break
      #} end if
    #} end for
    if (not extension_recognized): #{
      raise GroupFilteringError(
        "Unrecognized input extension: %s, " % input_file_name +
        "Valid Extensions: %s" % ", ".join(recognized_extensions))
    #} end if
    return os.path.join(self.options.output_dir, output_file_name)
  #} end def
#} end class

class FilterCls: #{
  def __init__(self, description, fail_reason, filter_type, value): #{
    self.description = description
    self.fail_reason = fail_reason
    self.filter_type  = filter_type
    self.value = value
  #} end def

  def ValueString(self): #{
    if (None == self.value): #{
      return "none"
    elif (FT_BOOL == self.filter_type):
      if (self.value): #{
        return "yes"
      else:
        return "no"
      #} end if
    elif (FT_INT == self.filter_type):
      return "%i" % self.value
    elif (FT_PERCENT == self.filter_type):
      return "%.1f" % self.value
    elif (FT_FLOAT == self.filter_type):
      return "%.3f" % self.value
    elif (FT_SET == self.filter_type):
      return ",".join(sorted(self.value))
    else:
      raise GroupFilteringError("unrecognized filter type: %s" %
        self.filter_type)
    #} end if
  #} end def
#} end class

def AdHocGetContigRepresented(member, log_info=None, robust=False):
  if (not member.gap): #{
    return member.TotalContigRepresented()
  #} end if
  # check that the meta-data could be parsed
  #if (None == member.meta_match):
  if (0 == len(member.meta_fields)): #{
    raise ACEventMemberError("could not parse gap meta data for "
      "member %s: %s" % (member.IDString(), member.meta_data))
  #} end if
  # ensure that the gap alignment does not extend outside of the gap
  gap_align_start = member.align_info_B.ctg_start
  gap_start       = member.meta_fields['gap_start']
  gap_align_start = max(gap_start, gap_align_start)
  gap_align_end   = member.align_info_B.ctg_end
  gap_end         = member.meta_fields['gap_end']
  gap_align_end   = min(gap_end, gap_align_end)
  aligned_gap_span = (gap_align_end - gap_align_start) + 1
  DebugMsg(log_info, "Contig align coords: %s-%s" %
    (gap_align_start, gap_align_end))
  # if the gap alignment did not fall inside the gap
  if (1 > aligned_gap_span and "dup" in member.topology): #{
    # try using the DUP coordinates
    gap_align_start = member.meta_fields['dup_start']
    gap_align_start = max(gap_start, gap_align_start)
    gap_align_end   = member.meta_fields['dup_end']
    gap_align_end   = min(gap_end, gap_align_end)
    aligned_gap_span = (gap_align_end - gap_align_start) + 1
    DebugMsg(log_info, "Invalid gap align span: %i, using DUP coords\n" %
      aligned_gap_span + "Gap align coords: %s-%s" %
      (gap_align_start, gap_align_end))
    #} end if
  #} end if
  # if the gap alignment still did not fall inside the gap
  if (1 > aligned_gap_span): #{
    msg_list = list([
      "could not get aligned gap span from member: %s" % member.IDString(),
      "CTG1:%s" % member.align_info_A.ContigCoordsString(),
      "CTG2:%s" % member.align_info_B.ContigCoordsString(),
      "DUP:%s"  % member.DupCoordsString(),
      "GAP:%s"  % member.GapCoordsString(),
      "aligned gap span: %i" % aligned_gap_span,
    ])
    if (robust): #{
      LogMsg(log_info, "WARNING: %s" % " ".join(msg_list))
    else:
      raise GroupFilteringError("\n  ".join(msg_list))
    #} end if
  #} end if
  # before_gap_span + after_gap_span + aligned_gap_span
  ctg_represented = (member.BeforeGapSpan() +
                     member.AfterGapSpan() +
                     aligned_gap_span)
  ctg_represented_fraction = (float(ctg_represented) /
                              float(member.contig_info.length))
  DebugMsg(log_info, "\n".join([
    "Before Gap: %i"    % member.BeforeGapSpan(),
    "After Gap: %i"     % member.AfterGapSpan(),
    "Gap: %i"           % aligned_gap_span,
    "Total: %i"         % ctg_represented,
    "Contig length: %i" % member.contig_info.length,
    "Fraction: %.3f"      % ctg_represented_fraction]))
  if (1 < ctg_represented_fraction): #{
    raise GroupFilteringError("Invalid contig representation fraction: "
      "%.3f for member:\n%s" % (ctg_represented_fraction, member.DataString()))
  #} end if
  return ctg_represented_fraction
#} end def

#def GetGeneNamesFromFile(path, log_info=None): #{
#  if (None == path): #{
#    return None
#  #} end if
#  start_time = time.time()
#  LogMsg(log_info, "Getting reference gene names...")
#  annots_type = GetAnnotationsType(path)
#  parse_annot = GetAnnotationsParseFunction(annots_type, log_info)
#  gene_name_set = set()
#  gene_names_file = FileBoxCls(path, "r", "cannot open gene names file")
#  for gene_line in gene_names_file: #{
#    annot = FixAnnotation(parse_annot(CleanLine(gene_line)))
#    # skip blank lines
#    if (None == annot): #{
#      continue
#    #} end if
#    try: #{
#      gene_name = annot.alias
#    except IndexError, e:
#      raise GroupFilteringError("error parsing line while getting "
#        "reference gene names:\n  %s" % CleanLine(gene_line) + "\n%s" % e)
#    #} end try
#    DebugMsg(log_info, "Gene name: %s" % gene_name)
#    #} end if
#    gene_name_set.update([gene_name.lower()])
#  #} end for
#  gene_names = sorted(gene_name_set)
#  LogMsg(log_info, "Time spent getting gene names: %s" %
#    TimeSpent(start_time))
#  return gene_names
#} end def

#### EXCEPTION CLASSES ####
class GroupFilteringError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Applies the given filters to the candidate contig "
    "groups that Barnacle identified.")
  args = [ "LIB", "BARNACLE_FILE", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--read-to-contig-with-weak",
                    type="int", dest="read_to_ctg", metavar="N",
                    help="Fail groups with no contig with more than N "
                         "supporting read-to-contig alignments, unless they "
                         "pass the strong-only read-to-contig filter. [default: "
                         "%default]")
  parser.add_option("--read-to-contig",
                    type="int", dest="read_to_ctg_u", metavar="N",
                    help="Fail groups with no contig with more than N "
                         "strongly-supporting read-to-contig alignments, unless "
                         "they pass the with-weak read-to-contig filter. "
                         "[default: %default]")
  parser.add_option("--read-to-contig-weak",
                    action="store_true", dest="req_u_read_to_ctg",
                    help="Do not fail groups failing the unique "
                         "read-to-contig filter, unless they also fail the "
                         "total read-to-contig filter. [default]")
  parser.add_option("--read-to-contig-strong",
                    action="store_true", dest="req_u_read_to_ctg",
                    help="Fail groups failing the unique read-to-contig "
                         "filter, even if they would pass the total "
                         "read-to-contig filter.")
  parser.add_option("--pair-to-genome",
                    type="int", dest="pair_to_gen", metavar="N",
                    help="Fail groups with fewer than N supporting "
                         "read-pair-to-genome alignments. [default: %default]")
  parser.add_option("--filtered-pair-to-genome",
                    type="int", dest="pair_to_gen_f", metavar="N",
                    help="Fail groups with fewer than N mapq-filtered "
                         "supporting read-pair-to-genome alignments. "
                         "[default: %default]")
  parser.add_option("--ignore-intronic",
                    action="store_true", dest="ignore_intron",
                    help="Do not count read-pairs mapping to inferred "
                         "intronic locations when calculating "
                         "read-pair-to-genome support. [default]")
  parser.add_option("--count-intronic",
                    action="store_false", dest="ignore_intron",
                    help="Count read-pairs mapping to inferred intronic "
                         "locations when calculating read-pair-to-genome "
                         "support.")
  parser.add_option("--no-struct-RNA",
                    action="store_true", dest="no_rna",
                    help="Fail groups overlapping small structural RNA "
                         "regions. [default]")
  parser.add_option("--allow-struct-RNA",
                    action="store_false", dest="no_rna",
                    help="Do not fail groups overlapping small structural "
                         "RNA regions.")
  parser.add_option("--no-mitochondria",
                    action="store_true", dest="no_mito",
                    help="Fail groups involving mitochondrial DNA. "
                         "[default]")
  parser.add_option("--allow-mitochondria",
                    action="store_false", dest="no_mito",
                    help="Do not fail groups involving mitochondrial DNA.")
  parser.add_option("--max-num-groups",
                    type="int", dest="max_num_groups", metavar="N",
                    help="Fail groups with no contig appearing in fewer than "
                         "N groups. Set this to zero to disable this filter. "
                         "[default: %default]")
  parser.add_option("--allow-single-aligner",
                    action="store_false", dest="reguire_multiple_aligns",
                    help="Do not fail groups based on the number of aligners "
                         "they were found with. [default]")
  parser.add_option("--require-multiple-aligners",
                    action="store_true", dest="reguire_multiple_aligns",
                    help="Fail groups that were not found with multiple "
                         "aligners.")
  parser.add_option("--min-identity",
                    type="float", dest="min_pid", metavar="F",
                    help="Fail groups with no contig that has greater than "
                         "F% percent identity for both alignments. [default: "
                         "%default]")
  parser.add_option("--length-sensitive-pid",
                    action="store_true", dest="adjust_pid",
                    help="Adjust the minimum percent identity to use "
                         "maximum number of mismatches for short alignment "
                         "blocks. [default]")
  parser.add_option("--absolute-pid",
                    action="store_false", dest="adjust_pid",
                    help="Use the minimum percent identity provided without "
                         "adjustment from maximum number of mismatches for "
                         "short alignment blocks.")
  parser.add_option("--max-mismatches",
                    type="int", metavar="N",
                    help="If --length-sensitive-pid is used, use the "
                         "minimum of the given minimum percent identity "
                         "value and the calculated percent identity of an "
                         "alignment with N mismatches. [default=%default]")
  parser.add_option("--no-homopolymers",
                    action="store_true", dest="no_runs",
                    help="Fail gapped events involving runs of a single "
                         "base. [default]")
  parser.add_option("--allow-homopolymers",
                    action="store_false", dest="no_runs",
                    help="Do not fail gapped events involving runs of a "
                         "single base.")
  parser.add_option("--soft-runs",
                    action="store_true", dest="soft_runs",
                    help="Allow a single divergent base when determining "
                         "whether a sequence will be considered a homopolymer "
                         "run. [default]")
  parser.add_option("--hard-runs",
                    action="store_false", dest="soft_runs",
                    help="Strictly require all bases to be identical for a "
                         "sequence to be considered a homopolymer run.")
  parser.add_option("--max-ctg-overlap",
                    type="int", dest="ctg_olap", metavar="N",
                    help="Fail paired-alignment candidates when the "
                         "contig coordinates of the alignments overlap by "
                         "more than Nbp. It makes sense for this value to be "
                         "approximately one read-length. [default: %default]")
  parser.add_option("--min-ctg-rep",
                    type="float", dest="ctg_rep", metavar="F",
                    help="Fail groups with no contig with the fraction of "
                         "its length represented by the alignments being at "
                         "least F. [default: %default]")
  parser.add_option("--no-multi-mapping",
                    action="store_true", dest="filter_multi_maps",
                    help="Fail groups with alignments that are flagged "
                         "as multi-mapping. [default]")
  parser.add_option("--allow-multi-mapping",
                    action="store_false", dest="filter_multi_maps",
                    help="Do not fail groups with alignments that are "
                         "flagged as multi-mapping.")
  parser.add_option("--no-repeats",
                    action="store_true", dest="filter_repeats",
                    help="Fail groups that are flagged as "
                         "overlapping repeat sequence.")
  parser.add_option("--allow-repeats",
                    action="store_false", dest="filter_repeats",
                    help="Do not fail groups that are flagged as "
                         "overlapping repeat sequence. [default]")
  parser.add_option("--no-polyA-events",
                    action="store_true", dest="filter_polyA",
                    help="Fail groups that are probably polyA "
                         "tails: single-base runs at the very "
                         "beginning or end of a contig. [default]")
  parser.add_option("--allow-polyA-events",
                    action="store_false", dest="filter_polyA",
                    help="Do not fail groups that are probably "
                         "polyA tails: single-base runs at the very "
                         "beginning or end of a contig.")
  parser.add_option("--ignore-topologies",
                    metavar="TOPOLOGIES",
                    help="Ignore groups classified as any of the topologies "
                         "in the comma-separated list TOPOLOGIES.")
  parser.add_option("--sort-by-r2cu",
                    action="store_true",
                    help="Output the groups sorted by their unique "
                         "read-to-contig coverage in descending order.")
  parser.add_option("-p", "--pretty",
                    action="store_true",
                    help="Print output in a more readable, but less easily "
                         "searched/parsed format.")
  parser.add_option("--data-check",
                    action="store_true",
                    help="Perform some simple sanity checks on the groups to "
                         "ensure that the upstream code is not having any "
                         "obvious problems. [default]")
  parser.add_option("--no-data-check",
                    action="store_false", dest="data_check",
                    help="Presume that the input file is well-formatted.")
  parser.add_option("--no-split-out",
                    action="store_false", dest="split_out",
                    help="Do not output groups of different topologies to "
                         "separate output files (create only pass and fail "
                         "output files). [default]")
  parser.add_option("--split-out",
                    action="store_true", dest="split_out",
                    help="Output groups with different alignment topologies "
                         "to separate output files.")
  parser.add_option("--print-fails",
                    action="store_true", dest="print_fails",
                    help="Write failing groups to a failed groups file. "
                         "[default]")
  parser.add_option("--no-print-fails",
                    action="store_false", dest="print_fails",
                    help="Do not write failing groups to a failed groups "
                         "file.")
  parser.add_option("--filter-gene-names",
                    dest="gene_names_path", metavar="FILE",
                    help="Use only gene names found in FILE (file must have "
                         "one gene name per line).")
  parser.add_option("--recalculate-ctg-rep",
                    action="store_true",
                    help="Recalculate the contig representation fraction "
                         "for gap events. [default]")
  parser.add_option("--no-recalculate-ctg-rep",
                    action="store_true",
                    help="Always use the reported contig representation "
                         "fraction.")
  parser.add_option("--robust-crc",
                    action="store_true", dest="robust_crc",
                    help="If a reasonable contig representation cannot be "
                         "calculated for a gap event, just fail that "
                         "candidate and display a warning rather than "
                         "raising an error. [default]")
  parser.add_option("--brittle-crc",
                    action="store_false", dest="robust_crc",
                    help="If a reasonable contig representation cannot be "
                         "calculated for a gap event, raise an error.")
  parser.add_option("--member-wise",
                    action="store_true", dest="member_wise",
                    help="At least one member of a group must pass all "
                         "filters for the whole group to pass. Only output "
                         "group members that pass all filters. [default]")
  parser.add_option("--group-wise",
                    action="store_false", dest="member_wise",
                    help="If any member of a group passes a filter, then the "
                         "whole group passes that filter.")
  misc_group = OptionGroup(parser, "Miscellaneous Options")
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
  parser.set_defaults(read_to_ctg=100,
                      read_to_ctg_u=5,
                      req_u_read_to_ctg=False,
                      pair_to_gen=0,
                      pair_to_gen_f=0,
                      ignore_intron=True,
                      no_rna=True,
                      no_mito=True,
                      max_num_groups=3,
                      reguire_multiple_aligns=False,
                      min_pid=99.0,
                      adjust_pid=True,
                      max_mismatches=1,
                      no_runs=True,
                      soft_runs=True,
                      ctg_olap=75,
                      ctg_rep=0.90,
                      filter_multi_maps=True,
                      filter_repeats=False,
                      filter_polyA=True,
                      ignore_topologies="gap-nontandem-inverted_duplication,gap-tandem-inverted_duplication,local-inversion",
                      sort_by_r2cu=False,
                      pretty=False,
                      check_data=True,
                      split_out=False,
                      print_fails=True,
                      recalculate_ctg_rep=True,
                      robust_crc=True,
                      member_wise=True,
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
  # get the input directory
  (input_dir, file_name) = os.path.split(options.barnacle_path)
  # get and check the output path
  options.output_dir = GetOutDir(input_dir, "filtered")
  # make sure we are not overwriting anything
  #if (not options.force and os.path.exists(options.output_dir)): #{
  #  path_errors.append("the output director already exists, please move "
  #    "the old directory so that it will not be overwritten or use the "
  #    "--force option to force overwriting.")
  #} end if
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "filter", options.output_dir)
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
    if (CheckPaths(options)): #{
      try: #{
        quick_filter = GroupFilteringCls(options)
        WriteCommand(quick_filter, sys.argv)
        quick_filter.FilterFile()
        quick_filter.DisplayFilterResults()
      except (MyError), e:
        ErrMsg("ERROR while running quick_filter:\n  %s" % e)
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
