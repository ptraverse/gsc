#! /usr/bin/env python
"""
candidate_group.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# TODO

# import standard modules
import re, string

# import custom modules
from utils.error import MyError
from utils.general import NormalizeChrID
from utils.messages import ErrMsg, DebugMsg
from parsers.tokenizer import TokenizerCls, ParseWarningCls, GetFieldAndValue
from grouped_candidate import GroupedCandidateCls, GroupedCandidateError
from coord_pair import CoordPairCls

# CONSTANTS

class CandidateGroupCls: #{
  def __init__(self, data_str): #{
    self.SetupParseFunctions()
    self.warnings = list()
    self.member_warnings = list()
    self.fail_reasons = set()
    self.members = list()
    self.candidates = self.members
    self.gap = False
    self.max_read_to_ctg_u = 0
    self.max_member_score = 0.0
    self.score = None
    self.total_members = None
    self.gene_names_filtered = False
    self.gene_dirs = None
    self.extra_data = list()
    self.rna = None
    try:
      self.ParseDataString(data_str)
    except GroupedCandidateError, e:
      raise CandidateGroupError("Error parsing data string for group:\n"
        "%s\n%s" % (data_str, e))
    # end try
  #} end def

  def SetupParseFunctions(self): #{
    self.ParseField = dict({
      'GRPNUM':         self.ParseGroupNumber,
      # TYPE is temporary, just until all files have it changed to EVENT_TYPE
      'TYPE':           self.ParseEventType,
      'EVENT_TYPE':     self.ParseEventType,
      'TOPOLOGIES':     self.ParseTopologies,
      'COORDS':         self.ParseCoords,
      'MEMBERS':        self.ParseNumMembers,
      'CONTIGS':        self.ParseNumContigs,
      'ALIGNERS':       self.ParseNumAligners,
      'CLENS':          self.ParseCtgLengths,
      'PAIR_TO_GENOME': self.ParsePairToGenome,
      #'EXON_BOUNDS':    self.ParseExonBounds,
      'RNA':            self.ParseRNA,
      'STATUS':         self.ParseStatus,
      'LIB':            self.ParseLib,
    })
  #} end def

  def __lt__(self, other): #{
    return self.score < other.score
  #} end def

  def ParseDataString(self, data_str): #{
    # tokenize the data string
    tokens = TokenizerCls(data_str, " ")
    for token in tokens: #{
      try:
        (field, value) = GetFieldAndValue(token)
        #ErrMsg("Parsing %s: %s" % (field, value))
      except ValueError, e:
        self.AddWarning("cannot separate field and value for token", token)
        continue
      # end try
      # if the field is a pair-to-genome field
      if (field.startswith("PAIR_TO_GENOME")): #{
        # move the mapq filter value from the field label to the field value
        (field, min_mapq) = field.rsplit("_", 1)
        value = (min_mapq, value)
        #ErrMsg("PAIR_TO_GENOME! Parsing %s: %s" % (field, value))
      #} end if
      # if the field is a standard field
      if (field in self.ParseField): #{
        self.ParseField[field](value)
      else:
        #ErrMsg("Parsing extra field: %s" % field)
        # it is just an extra field
        self.extra_data.append(token)
      #} end if
    #} end for
    if (0 < len(self.extra_data)): #{
      self.tail = " ".join(self.extra_data)
    #} end if
    #self.ParseExtraInfo(extra_info)
  #} end def

  def ParseGroupNumber(self, group_num_str): #{
    try:
      self.id = int(group_num_str)
    except ValueError, e:
      self.AddWarning("invalid group ID", group_num_str, "must be an integer")
    # end try
  #} end def

  def ParseEventType(self, event_type_str): #{
    self.event_type = event_type_str
  #} end def

  def ParseTopologies(self, topologies_str): #{
    self.topologies = sorted(set(topologies_str.split(",")))
  #} end def

  def ParseCoords(self, coords_str): #{
    # for split candidates
    if ";" in coords_str: #{
      (coord_str_A, coord_str_B) = coords_str.split(";")
      self.coord_pair_A = ChromAndCoordPairCls(coord_str_A)
      self.coord_pair_B = ChromAndCoordPairCls(coord_str_B)
    # for gap candidates
    else:
      self.coord_pair_A = ChromAndCoordPairCls(coords_str)
      self.coord_pair_B = None
    #} end if
  #} end def

  def ParseNumMembers(self, num_members_str): #{
    try:
      if ("(" in num_members_str): #{
        num_members_patt = r"\A(?P<filtered>\d+)\((?P<total>\d+)\)\Z"
        num_members_match = re.search(num_members_patt, num_members_str)
        if (None == num_members_match): #{
          self.AddWarning("unable to parse number of members", num_members_str)
          return
        #} end if
        self.num_members = int(num_members_match.group('filtered'))
        self.total_members = int(num_members_match.group('total'))
      else:
        self.num_members = int(num_members_str)
      #} end if
    except ValueError, e:
      self.AddWarning("invalid number of members", num_members_str,
        "must be an integer")
    # end try
  #} end def

  def ParseNumContigs(self, num_contigs_str): #{
    try:
      self.num_contigs = int(num_contigs_str)
    except ValueError, e:
      self.AddWarning("invalid number of contigs", num_contigs_str,
        "must be an integer")
    # end try
  #} end def

  def ParseNumAligners(self, num_aligners_str): #{
    try:
      self.num_aligners = int(num_aligners_str)
    except ValueError, e:
      self.AddWarning("invalid number of aligners", num_aligners_str,
        "must be an integer")
    # end try
  #} end def

  def ParseCtgLengths(self, ctg_lengths_str): #{
    try:
      self.ctg_lengths = [int(length.replace("bp","")) for
        length in ctg_lengths_str.split(",")]
    except ValueError, e:
      self.AddWarning("invalid contig length", ctg_lengths_str,
        "must be integers")
    # end try
  #} end def

  def ParsePairToGenome(self, pair_to_genome_data): #{
    (min_mapq_str, pair_to_genome_str) = pair_to_genome_data
    #ErrMsg("Min MapQ: " + min_mapq_str + " P2G: " + pair_to_genome_str)
    try:
      min_mapq = int(min_mapq_str)
    except ValueError, e:
      self.AddWarning("invalid minimum mapq value", min_mapq_str,
        "must be an integer")
      return
    # end try
    if (0 == min_mapq): #{
      self.ParsePairToGenomeTotal(pair_to_genome_str)
    else:
      self.ParsePairToGenomeFiltered(min_mapq, pair_to_genome_str)
    #} end if
  #} end def

  def ParsePairToGenomeTotal(self, pair_to_genome_str): #{
    (self.pair_to_genome_exonic,
     self.pair_to_genome_all) = self.GetPairToGenomeValues(pair_to_genome_str)
    if (self.pair_to_genome_exonic > self.pair_to_genome_all): #{
      self.AddWarning("invalid pair-to-genome values", pair_to_genome_str,
        "exonic support must be less than total pair support")
    #} end if
  #} end def

  def ParsePairToGenomeFiltered(self, min_mapq,
      pair_to_genome_filtered_str): #{
    self.min_mapq = min_mapq
    (self.pair_to_genome_filtered_exonic,
     self.pair_to_genome_filtered_all) = self.GetPairToGenomeValues(
        pair_to_genome_filtered_str)
    if (self.pair_to_genome_filtered_exonic >
        self.pair_to_genome_filtered_all): #{
      self.AddWarning("invalid pair-to-genome values",
        pair_to_genome_filtered_str,
        "exonic support must be less than total pair support")
    #} end if
  #} end def

  def GetPairToGenomeValues(self, pair_to_genome_str): #{
    if ("(" in pair_to_genome_str): #{
      p2g_pattern = r"(?P<exon>\d+)\((?P<all>\d+)\)"
      p2g_match = re.search(p2g_pattern, pair_to_genome_str)
      if (None == p2g_match): #{
        self.AddWarning("cannot parse pair to genome value", pair_to_genome_str)
        return (None, None)
      #} end if
      exonic = int(p2g_match.group('exon'))
      all    = int(p2g_match.group('all'))
      return (exonic, all)
    elif ("N/A" == pair_to_genome_str): #{
      return (pair_to_genome_str, pair_to_genome_str)
    else:
      try:
        pair_to_genome = int(pair_to_genome_str)
      except (TypeError, ValueError), e:
        self.AddWarning("invalid pair to genome value",
          pair_to_genome_str, "must be an integer")
        return (None, None)
      # end try
      return (pair_to_genome, pair_to_genome)
    #} end if
  #} end def

  def ParseExonBounds(self, exon_bounds_str): #{
    if ("N/A" == exon_bounds_str): #{
      self.exon_bounds = None
    elif (exon_bounds_str in ["0", "1", "2"]): #{
      self.exon_bounds = int(exon_bounds_str)
    else:
      self.AddWarning("invalid exon bounds count", exon_bounds_str,
        "must be 0, 1, or 2")
    #} end if
  #} end def

  def ParseRNA(self, rna_str): #{
    if ("N/A" == rna_str): #{
      self.rna = None
    elif (rna_str in ["Y", "N"]): #{
      self.rna = ("Y" == rna_str)
    else:
      self.AddWarning("invalid RNA overlap value", rna_str,
        "must be \"Y\" or \"N\"")
      return
    #} end if
  #} end def

  def ParseStatus(self, status_str): #{
    if ("RNA_REPEAT" in status_str): #{
      self.rna = True
    #} end if
    if (status_str.startswith("FAIL")): #{
      fail_pattern = r"FAIL\((?P<fail_reasons>.+)\)"
      fail_match = re.search(fail_pattern, status_str)
      if (None == fail_match): #{
        self.AddWarning("cannot parse fail reasons", status_str)
        return
      #} end if
      self.fail_reasons = set(fail_match.group('fail_reasons').split(","))
    #} end if
  #} end def

  def ParseLib(self, lib_str): #{
    self.lib = lib_str
  #} end def

  def DataList(self, readable=False): #{
    data_list = [
      "GRPNUM:%i" % self.id,
    ]
    if (hasattr(self, "event_type")): #{
      data_list.append("EVENT_TYPE:%s" % self.event_type)
    #} end if
    data_list.extend([
      "TOPOLOGIES:%s" % ",".join(self.topologies),
      "COORDS:%s" % self.CoordsString(),
      "MEMBERS:%s" % self.NumMembersString(readable),
      "CONTIGS:%i" % self.num_contigs,
      #"ALIGNERS:%i" % self.num_aligners,
      "CLENS:%s" % self.CtgLengthsString(),
      "PAIR_TO_GENOME_0:%s" % self.PairToGenome(),
      self.PairToGenomeFiltered(),
      #"EXON_BOUNDS:%s" % self.ExonBoundsString(),
      "RNA:%s" % self.RNAString()
    ])
    data_list.append("STATUS:%s" % self.Status())
    if (hasattr(self, "lib") and None != self.lib): #{
      data_list.append("LIB:%s" % self.lib)
    #} end if
    data_list.extend(self.extra_data)
    if (None != self.score): #{
      data_list.append("SCORE:%.2f" % self.score)
    #} end if
    return data_list
  #} end def

  def DataString(self): #{
    data_fields = self.DataList()
    return " ".join(data_fields)
  #} end def

  def ReadableString(self): #{
    data_fields = ["="*80]
    data_fields.extend(self.DataList(readable=True))
    if (None != self.score): #{
      # insert the group score field
      data_fields.insert(2, data_fields[-1])
      del data_fields[-1]
    #} end if
    return "\n".join(data_fields)
  #} end def

  def CoordsString(self): #{
    if (None == self.coord_pair_B): #{
      return self.coord_pair_A.ToString()
    else:
      coord_str_A = self.coord_pair_A.ToString()
      coord_str_B = self.coord_pair_B.ToString()
      return ";".join([coord_str_A, coord_str_B])
    #} end if
  #} end def

  def NumMembersString(self, readable=False): #{
    num_member_str = "%i" % self.num_members
    if (None != self.total_members): #{
      if (readable): #{
        num_member_str += " TOTAL_MEMBERS:%i" % self.total_members
      else:
        num_member_str += "(%i)" % self.total_members
      #} end if
    #} end if
    return num_member_str
  #} end def

  def CtgLengthsString(self): #{
    return ",".join(["%ibp" % length for length in self.ctg_lengths])
  #} end def

  def PairToGenome(self): #{
    if ("N/A" == self.pair_to_genome_exonic): #{
      return self.pair_to_genome_exonic
    #} end if
    if (self.pair_to_genome_exonic < self.pair_to_genome_all): #{
      return "%i(%i)" % (self.pair_to_genome_exonic, self.pair_to_genome_all)
    else:
      return "%i" % (self.pair_to_genome_exonic)
    #} end if
  #} end def

  def PairToGenomeFiltered(self): #{
    if ("N/A" == self.pair_to_genome_filtered_exonic): #{
      p2g_str = self.pair_to_genome_filtered_exonic
    elif (self.pair_to_genome_filtered_exonic <
          self.pair_to_genome_filtered_all): #{
      p2g_str = "%i(%i)" % (self.pair_to_genome_filtered_exonic,
        self.pair_to_genome_filtered_all)
    else:
      p2g_str = "%i" % (self.pair_to_genome_filtered_exonic)
    #} end if
    return "PAIR_TO_GENOME_%i:%s" % (self.min_mapq, p2g_str)
  #} end def

  def ExonBoundsString(self): #{
    if (None == self.exon_bounds): #{
      return "N/A"
    else:
      return "%i" % self.exon_bounds
    #} end if
  #} end def

  def RNAString(self): #{
    if (None == self.rna): #{
      return "N/A"
    #} end if
    elif (self.rna): #{
      return "Y"
    else:
      return "N"
    #} end if
  #} end def

  def Status(self): #{
    if (0 < len(self.fail_reasons)): #{
      return "FAIL(%s)" % ",".join(sorted(self.fail_reasons))
    else:
      return "PASS"
    #} end if
  #} end def

  def PassedFilters(self): #{
    if (0 < len(self.fail_reasons)): #{
      return False
    else:
      return True
    #} end if
  #} end def

  def AddWarning(self, warning, data, rule=""): #{
    parse_warning = ParseWarningCls(warning, data, rule)
    self.warnings.append(parse_warning)
  #} end def

  def WarningsString(self, indent=""): #{
    warnings_list = list()
    for parse_warning in self.warnings: #{
      message = ("Warning: %s for group %i: \"%s\"" %
        (parse_warning.warning, self.id, parse_warning.data))
      if ("" != parse_warning.rule): #{
        message += ", %s" % parse_warning.rule
      #} end if
      warnings_list.append(message)
    #} end for
    delim = "%s\n" % indent
    return delim.join(warnings_list)
  #} end def

  def AddMemberFromString(self, member_data_str, check_data=False): #{
    try:
      new_member = GroupedCandidateCls(member_data_str, check_data)
    except GroupedCandidateError, e:
      raise CandidateGroupError("Error adding member to group:\n%s\n%s\n%s" %
        (self.DataString(), member_data_str, e))
    # end try
    if (new_member.group_id != self.id): #{
      raise CandidateGroupError("Member group id \"%i\" " % member.group_id +
        "does not equal group id \"%i\"." % group.id)
    #} end if
    self.members.append(new_member)
    if (0 < len(new_member.warnings)): #{
      self.member_warnings.append(new_member.WarningsString())
    #} end if
    if (new_member.gap): #{
      self.gap = True
    #} end if
    if (new_member.read_to_ctg_unique > self.max_read_to_ctg_u): #{
      self.max_read_to_ctg_u = new_member.read_to_ctg_unique
    #} end if
  #} end def

  def UpdateMemberCount(self, new_num_members=None): #{
    if (None == self.total_members): #{
      self.total_members = self.num_members
    #} end if
    if (None == new_num_members): #{
      new_num_members = len(self.members)
    #} end if
    self.num_members = new_num_members
  #} end def

  def RemoveLowScoringMembers(self, min_member_score): #{
    new_members = list()
    contigs = set()
    contig_lens = set()
    aligners = set()
    for member in self.members: #{
      if (member.score >= min_member_score): #{
        new_members.append(member)
        contigs.add(member.contig_info.id)
        if (member.aligner.startswith("gap")): #{
          aligners.add("blat")
        else:
          aligners.add(member.aligner.split(",", 1)[0])
        #} end if
        contig_lens.add(member.contig_info.length)
      #} end if
    #} end for
    self.members = new_members
    self.UpdateMemberCount()
    self.num_contigs  = len(contigs)
    self.num_aligners = len(aligners)
    self.ctg_lengths  = sorted(contig_lens, reverse=True)
  #} end def

  def FullDataString(self, readable=False): #{
    if (readable): #{
      full_data_list = [self.ReadableString()]
    else:
      full_data_list = [self.DataString()]
    #} end if
    for member in self.members: #{
      if (readable): #{
        full_data_list.append(member.ReadableString())
      else:
        full_data_list.append(member.DataString())
      #} end if
    #} end for
    return "\n".join(full_data_list)+"\n\n"
  #} end def

  def ClearBPGenes(self): #{
    for member in self.members: #{
      member.ClearBPGenes()
    #} end for
  #} end def

  def GeneSet(self, retain_case=False): #{
    gene_set = set()
    for member in self.members: #{
      if (member.gap): #{
        member_list = member.genes_A
      else:
        member_list = member.genes_A + member.genes_B
      #} end if
      if (retain_case): #{
        gene_set.update(member_list)
      else:
        gene_set.update([gene.lower() for gene in member_list])
      #} end if
    #} end for
    return gene_set
  #} end def

  def GetGeneSets(self): #{
    genesA = set(self.members[0].genes_A)
    if (self.members[0].gap): #{
      genesB = set()
    else:
      genesB = set(self.members[0].genes_B)
    #} end if
    if (1 == len(self.members)): #{
      return (genesA, genesB)
    #} end if
    for member in self.members[1:]: #{
      overlap = (len(genesA.intersection(member.genes_A)) +
        len(genesB.intersection(member.genes_B)))
      swap_overlap = (len(genesA.intersection(member.genes_B)) +
        len(genesB.intersection(member.genes_A)))
      if (overlap >= swap_overlap): #{
        genesA.update(member.genes_A)
        genesB.update(member.genes_B)
      else:
        genesA.update(member.genes_B)
        genesB.update(member.genes_A)
      #} end if
    #} end for
    return (genesA, genesB)
  #} end def

  def FilterGeneNames(self, ref_gene_names): #{
    if (None == ref_gene_names or self.gene_names_filtered): #{
      return
    #} end if
    good_genes = self.GeneSet().intersection(ref_gene_names)
    for member in self.members: #{
      filtered_genes_A = list()
      for gene in member.genes_A: #{
        if (gene.lower() in good_genes): #{
          filtered_genes_A.append(gene)
        #} end if
      #} end for
      if (0 < len(filtered_genes_A)): #{
        member.genes_A = filtered_genes_A
      #} end if
      filtered_genes_B = list()
      for gene in member.genes_B: #{
        if (gene.lower() in good_genes): #{
          filtered_genes_B.append(gene)
        #} end if
      #} end for
      if (0 < len(filtered_genes_B)): #{
        member.genes_B = filtered_genes_B
      #} end if
    #} end for
    self.gene_names_filtered = True
  #} end def

  def TopologyString(self): #{
    if (1 == len(self.topologies)): #{
      return self.topologies[0]
    #} end if
    return 'multi'
  #} end def

  def ClearRepeats(self): #{
    for member in self.members: #{
      member.ClearRepeats()
    #} end for
  #} end def

  def AddRepeatsToGroup(self, repeats_dict, keep_current=False): #{
    for member_id in repeats_dict.keys(): #{
      #ErrMsg("Adding repeats to member: %s" % member_id)
      member_index = string.ascii_lowercase.find(member_id)
      #ErrMsg("  Index: %i" % member_index)
      repeats_dict_for_member = repeats_dict[member_id]
      self.members[member_index].AddRepeats(
        repeats_dict_for_member, keep_current)
    #} end for
  #} end def

  def OverlapsRepeat(self): #{
    for member in self.members: #{
      if (member.OverlapsRepeat()): #{
        return True
      #} end if
    #} end for
    return False
  #} end def

  def UpdateTopologies(self): #{
    topology_set = set()
    for member in self.members: #{
      if (member.PassedFilters()): #{
        topology_set.add(member.topology)
      #} end if
    #} end for
    self.topologies = sorted(topology_set)
  #} end def
#} end class

class ChromAndCoordPairCls(CoordPairCls): #{
  def __init__(self, coord_pair_str): #{
    self.ParseCoordPair(coord_pair_str)
  #} end def

  #def ParseCoordPair(self, coord_pair_str): #{}
  def ParseCoordPair(self, full_coord_string): #{
    #coord_pattern = r"chr(?P<chrom>\d+|[XY]|MT?):(?P<left>\d+)-(?P<right>\d+)"
    #coord_pattern = (r"(?P<chrom>%s):(?P<left>\d+)-(?P<right>\d+)" %
    #  CHR_ID_PATT)
    #coord_match   = re.search(coord_pattern, coord_pair_str)
    #if (None == coord_match): #{
    #  raise CandidateGroupError("cannot parse coordinate pair pattern from "
    #    "\"%s\"" % coord_pair_str)
    #} end if
    #self.chrom    = NormalizeChrID(coord_match.group('chrom'))
    #self.left     = int(coord_match.group('left'))
    #self.right    = int(coord_match.group('right'))
    (raw_chrom, coord_string) = full_coord_string.split(":", 1)
    self.chrom    = NormalizeChrID(raw_chrom)
    CoordPairCls.__init__(self, coord_string)
    (self.left, self.right) = (self.start, self.end)
    # the event's genomic coordinates should be ordered left-to-right
    if (self.right < self.left): #{
      raise CandidateGroupError("incorrectly ordered genomic coordinates")
    #} end if
  #} end def

  # Cannot determine strand from these coordinates,
  #  as they should always be ordered left-to-right
  #def Strand(self): #{
  #  if (self.start < self.end): #{
  #    return "+"
  #  else:
  #    return "-"
  #  #} end if
  #} end def

  def ToString(self): #{
    return "chr%s:%i-%i" % (self.chrom, self.left, self.right)
  #} end def
#} end class

class GroupParserCls: #{
  def __init__(self, keep_lines=False, log_info=None): #{
    self.keep_lines = keep_lines
    self.log_info = log_info
  #} end def

  def ParseGroup(self, group_line, group_data_file, check_data=False): #{
    #group_line = CleanLine(group_line)
    group_line = group_line
    DebugMsg(self, "Parsing group:\n%s" % group_line)
    group = CandidateGroupCls(group_line)
    if (self.keep_lines): #{
      self.group_line = group_line
      self.member_lines = list()
    #} end if
    for i in range(group.num_members): #{
      #member_line = CleanLine(group_data_file.next())
      member_line = group_data_file.next()
      DebugMsg(self, member_line)
      group.AddMemberFromString(member_line, check_data)
      if (self.keep_lines): #{
        self.member_lines.append(member_line)
      #} end if
    #} end def
    return group
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class CandidateGroupError(MyError): #{
  pass
#} end class
