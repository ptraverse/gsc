#! /usr/bin/env python
"""
fusion_predictor.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import re, time

# import custom modules
from utils.error import MyError
from utils.general import (TimeSpent, RemoveExonTypeLabel, ReverseComplement,
  NonStandardChr)
from utils.messages import LogMsg, DebugMsg, ExtremeDebugMsg
from utils.multi_dict import MultiDictCls
from utils.files_paths import GetFilePath, FileBoxCls
from base_predictor import BasePredictorCls
from parsers.genes.annotation import GeneAnnotationParserCls
from annotation.create_gene_feature_coords import INTRON

# constants
MIN_IMPROVE_FRACTION = 1.5

class FusionPredictorCls(BasePredictorCls): #{
  def __init__(self, options, log_info=None): #{
    BasePredictorCls.__init__(self, 'fusion', 'fus', 'fusion',
      options, log_info=log_info)
    #dup_gene_path = GetFilePath(options.output_dir, options.lib,
    #  "%s.dup_genes" % self.ext)
    #self.dup_gene_file = FileBoxCls(dup_gene_path, "w",
    #  "cannot create duplicated gene %s output file" % self.description)
    #mix_dir_path = GetFilePath(options.output_dir, options.lib,
    #  "%s.mix_dirs" % self.ext)
    #self.mix_dir_file = FileBoxCls(mix_dir_path, "w",
    #  "cannot create mixed gene directions %s output file" % self.description)
    self.AddOutFile("dup_gene", "%s.dup_genes" % self.ext,
      "duplicated gene %s" % self.description)
    self.AddOutFile("mix_dirs", "%s.mix_dirs" % self.ext,
      "mixed sense %s" % self.description)
    # TEMPORARY file to keep track of contigs for which one side is missing
    # from the realignment results
    self.AddOutFile("align_miss", "%s.align_miss.tmp" % self.ext,
      "missing realignment %s" % self.description)
    # file for potential "read-through" type events
    self.AddOutFile("readthru", "%s.rth" % self.ext,
      "read-through %s" % self.description)
    # gene_directions[id_type][id] = {+,-,*} where "*" = both directions
    self.gene_directions = {'gene_name':dict(), 'transcript_id':dict()}
    # get gene directions from annotation file
    if (self.options.use_gene_directions): #{
      self.GetGeneDirections()
    #} end if
    self.bad_topologies = set([
      'end-duplication',
      'junction-duplication',
    ])
  #} end def

  def GetGeneDirections(self): #{
    LogMsg(self, "Getting gene directions from annotations file...")
    start_time = time.time()
    DebugMsg(self, "  %s" % self.options.transcript_annotations)
    conflicts = {'gene_name':set(), 'transcript_id':set()}
    # get the appropriate parser
    #annot_type = GetAnnotationsType(self.options.transcript_annotations)
    #parse_annot = GetAnnotationsParseFunction(annot_type, self.log_info)
    # open the annotations file
    #annots_file = FileBoxCls(self.options.transcript_annotations, "r",
    #  "cannot open annotations file")
    # iterate through the annotations file
    #for annot_line in annots_file: #{}
    #  annot = parse_annot(annot_line)
    parser = GeneAnnotationParserCls(self.options.transcript_annotations,
      log_info=self.log_info)
    for annot in parser: #{
      # replace spaces in the name and alias
      #annot.name = annot.name.replace(" ", "_")
      #annot.alias = annot.alias.replace(" ", "_")
      #annot_name = annot.name.lower()
      #annot_alias = annot.alias.lower()
      # skip genes with already conflicting annotations
      #if (annot_name  in conflicted_genes or
      #    annot_alias in conflicted_genes): #{}
      if (annot.gene_name.lower()     in conflicts['gene_name'] and
          annot.transcript_id.lower() in conflicts['transcript_id']): #{
        continue
      #} end if
      # ensure that the gene is on a "normal" chromosome
      #chr_patt = r"^(chr)?([1-9][0-9]*|[XY]|MT?)$"
      #if (None == re.search(chr_patt, annot.chrom)): #{}
      if (NonStandardChr(annot.chrom)): #{
        continue
      #} end if
      #for field in ["name", "alias"]: #{}
      for field in ["transcript_id", "gene_name"]: #{
        id = getattr(annot, field).lower()
        if (id in self.gene_directions[field]): #{
          if (annot.strand != self.gene_directions[field][id]): #{
            msg = ("conflicting gene annotations: %s: %s, %s" %
              (getattr(annot, field), annot.strand,
              self.gene_directions[field][id]))
            # add the current gene to the list of conflicted genes
            #conflicted_genes.update([annot_name, annot_alias])
            conflicts[field].add(id)
            if (self.options.ignore_conflicts): #{
              # remove it from the directions dictionary
              del self.gene_directions[field][id]
              #if (annot_name in self.gene_directions): #{
              #  del self.gene_directions[annot_name]
              #} end if
              #if (annot_alias in self.gene_directions): #{
              #  del self.gene_directions[annot_alias]
              #} end if
              DebugMsg(self, "warning: %s" % msg)
            elif (self.options.use_conflicts):
              # mark the gene as going both ways
              self.gene_directions[field][id] = "*"
              #ExtremeDebugMsg(self, "warning: %s" % msg)
            else:
              raise FusionPredictorError(msg)
            #} end if
          #} end if
        else:
          self.gene_directions[field][id] = annot.strand
        #} end if
      #} end for
    #} end for
    #annots_file.Close()
    if (0 < len(conflicts['gene_name']) or
        0 < len(conflicts['transcript_id'])): #{
      LogMsg(self, "WARNING: %i genes (%i transcripts) are annotated as "
        "going both ways." % (len(conflicts['gene_name']),
        len(conflicts['transcript_id'])))
    #} end if
    LogMsg(self, "Time spent getting gene directions: %s" %
      TimeSpent(start_time))
  #} end def

  def TestGroup(self, group): #{
    DebugMsg(self, "Group %s. Topologies: %s" %
      (group.id, ",".join(group.topologies)))
    # gapped candidates are never fusions
    if (group.gap): #{
      DebugMsg(self, "gap candidates are never fusions")
      return False
    #} end if
    # interchromosomal candidates are always fusions
    # unless gene direction is being considered
    #if ("interchr" in group.topologies): #{
    #  if (0 == len(self.gene_directions)): #{
    #    DebugMsg(self, "interchromosomal candidate without considering "
    #      "gene directions")
    #    return True
    #  #} end if
    #else:
    # if the candidate is not interchromosomal, check whether
    # the genomic coordinates overlap
    if ("interchr" not in group.topologies): #{
      # if the candidate's genomic coordinates overlap,
      # then it is not a fusion
      if ((group.coord_pair_A.left <= group.coord_pair_B.left and
           group.coord_pair_B.left <  group.coord_pair_A.right) or
          (group.coord_pair_B.left <  group.coord_pair_A.left and
           group.coord_pair_A.left <  group.coord_pair_B.right)):
        DebugMsg(self, "genomic coordinates overlap")
        return False
      #} end if
      DebugMsg(self, "non-overlapping genomic coordinates")
    #} end if
    #return self.FusionGeneCheck(warnings)
    return True
  #} end def

  def TestMember(self, member): #{
    DebugMsg(self, "  Member %s. Topology: %s" %
      (member.candidate_id, member.topology))
    # do not call fusions with end or junction duplication topologies
    if (member.topology in self.bad_topologies): #{
      DebugMsg(self, "  %s: Wrong alignment topology: %s" %
        (member.IDString(), member.topology))
      return False
    #} end if
    # if either gene set is empty, then the event is not a fusion
    has_genes = dict()
    for region in ["A", "B"]: #{
      has_genes[region] = member.OverlapsGene(region,
        include_introns=self.options.include_introns)
      if (self.options.include_nearby and
          member.OverlapsGene("nearby_%s" % region)): #{
        has_genes[region] = True
      #} end if
    #} end for
    if (not has_genes["A"] or not has_genes["B"]): #{
      DebugMsg(self, "  empty gene set")
      return False
    #} end if
    # check the minimum unique alignment length threshold
    unique_A = (member.align_info_A.ContigSpan() -
      member.meta_fields['ctg_overlap'])
    unique_B = (member.align_info_B.ContigSpan() -
      member.meta_fields['ctg_overlap'])
    if (self.options.min_unique_align > unique_A or
        self.options.min_unique_align > unique_B): #{
      DebugMsg(self, "insufficient unique alignment length")
      return False
    #} end if
    # check the minimum exon boundaries threshold
    if (self.options.min_exon_bounds > member.meta_fields['exon_bounds']): #{
      DebugMsg(self, "insufficient exon boundaries")
      return False
    #} end if
    # if the candidate contig is not interchromosomal,
    # check that the gene sets are disjoint
    if ("interchr" != member.topology): #{
      DebugMsg(self, "  Genes: %s" % member.MainGenesString())
      # if the gene sets overlap at all, then the event is not a fusion
      if (member.GeneSetsOverlap()): #{
        DebugMsg(self, "  overlapping gene sets")
        return False
      #} end if
      DebugMsg(self, "  disjoint gene sets")
    #} end if
    # check gene directions and set 5' and 3' ends
    ExtremeDebugMsg(self, "Before gene direction check: %s" %
      member.MainGenesString())
    result = self.CheckGeneDirections(member)
    ExtremeDebugMsg(self, "After gene direction check: %s" %
      member.MainGenesString())
    return result
    #if (None != self.gene_directions): #{
    #  # check gene directions
    #  if (member.CheckGeneDirections(self.gene_dirs)): #{
    #  #} end if
    #} end if
  #} end def

  def CheckGeneDirections(self, member): #{
    DebugMsg(self, "  checking gene directions...")
    #if (None == self.gene_directions or
    #    not self.options.use_gene_directions): #{}
    if (0 == len(self.gene_directions['gene_name']) and
        0 == len(self.gene_directions['transcript_id'])): #{
      DebugMsg(self, "  No gene directions to check.")
      return True
    #} end if
    if (member.OverlapsGene("breakpoint_A", include_introns=True)): #{
      genes_A = HandleIntrons(member.breakpoint_genes_A, log_info=self.log_info)
    elif (member.OverlapsGene("A", include_introns=True)):
      genes_A = HandleIntrons(member.genes_A, log_info=self.log_info)
    else:
      genes_A = member.nearby_A
    #} end if
    orient_A = self.GetOrientation(genes_A,
      member.align_info_A.Strand())
    if ("?" == orient_A): #{
      if (member.OverlapsGene("A", include_introns=True) or
          member.OverlapsGene("nearby_A")):
        DebugMsg(self, "  could not get first gene direction for candidate.\n"
          "  %s" % member.DataString())
      #} end if
      return False
    #} end if
    if (member.OverlapsGene("breakpoint_B", include_introns=True)): #{
      genes_B = HandleIntrons(member.breakpoint_genes_B, log_info=self.log_info)
    elif (member.OverlapsGene("B", include_introns=True)):
      genes_B = HandleIntrons(member.genes_B, log_info=self.log_info)
    else:
      genes_B = member.nearby_B
    #} end if
    orient_B = self.GetOrientation(genes_B,
      member.align_info_B.Strand())
    if ("?" == orient_B): #{
      if (member.OverlapsGene("B", include_introns=True) or
          member.OverlapsGene("nearby_B")):
        DebugMsg(self, "  could not get second gene direction for candidate.\n"
          "  %s" % member.DataString())
      #} end if
      return False
    #} end if
    DebugMsg(self, "  Orientation A: %s, B: %s" % (orient_A, orient_B))
    if (("+" == orient_A and "+" == orient_B) or
        ("+" == orient_A and "*" == orient_B) or
        ("*" == orient_A and "+" == orient_B) or
        ("*" == orient_A and "*" == orient_B)): #{
      member.ends = ("5", "3")
    elif (("-" == orient_A and "-" == orient_B) or
        ("-" == orient_A and "*" == orient_B) or
        ("*" == orient_A and "-" == orient_B)):
      member.ends = ("3", "5")
    elif ("+" == orient_A and "-" == orient_B):
      member.ends = ("5", "5")
    elif ("-" == orient_A and "+" == orient_B):
      member.ends = ("3", "3")
    else:
      raise FusionPredictorError("cannot determine 5' and 3' ends of "
        "fusion, with orientations A: %s, B: %s" % (orient_A, orient_B))
    #} end if
    DebugMsg(self, "  Member ends: %s',%s'" % member.ends)
    if ("*" == orient_A or "*" == orient_B or orient_A == orient_B): #{
      DebugMsg(self, "  valid gene directions")
      return True
    #} end if
    DebugMsg(self, "non-continuous gene directions")
    # note: non-continuous gene direction (mixed-direction) events
    # will be stored in a separate file
    return True
    #return False
  #} end def

  def GetOrientation(self, gene_list, align_strand): #{
    # use the strand of the first gene found from the annotations file
    #plus_strand = False
    #minus_strand = False
    strand_flags = {'plus':False, 'minus':False}
    for gene in gene_list: #{
      # use transcripts if possible
      if (isinstance(gene_list, MultiDictCls) and
          2 < gene_list.depth): #{
        for transcript in gene_list[gene]: #{
          self.GetFeatureStrands(transcript, strand_flags, "transcript_id")
        #} end for
      else:
        self.GetFeatureStrands(gene, strand_flags, "gene_name")
      #} end if
    #} end for
    # if there are genes going both ways
    if (strand_flags['plus'] and strand_flags['minus']): #{
      return "*"
    # if the gene is on the same strand as the alignment
    elif ((strand_flags['plus']  and "+" == align_strand) or
          (strand_flags['minus'] and "-" == align_strand)):
      return "+"
    # if the gene is on the opposite strand as the alignment
    elif ((strand_flags['plus']  and "-" == align_strand) or
          (strand_flags['minus'] and "+" == align_strand)):
      return "-"
    # if the gene strand could not be determined (not strand_flags['plus'] and
    #   not strand_flags['minus']) or align_strand not in "+-"
    else:
      return "?"
    #} end if
    # use the strand of the first gene found from the annotations file
    #gene_strand = None
    #for gene in gene_list: #{
    #  gene_id = RemoveExonTypeLabel(gene).lower()
    #  if (gene_id in self.gene_directions): #{
    #    gene_strand = self.gene_directions[gene_id]
    #    break
    #  #} end if
    #} end for
    # if the gene strand could not be determined
    #if (None == gene_strand): #{
    #  return "?"
    # if the gene goes both ways
    #elif ("*" == gene_strand):
    #  return "*"
    # if the gene is on the same strand as the alignment
    #elif (gene_strand == align_strand):
    #  return "+"
    # if the gene is on the opposite strand as the alignment
    #else:
    #  return "-"
    #} end if
  #} end def

  def GetFeatureStrands(self, feature, strand_flags, field): #{
    feature_id = RemoveExonTypeLabel(feature).lower()
    if (feature_id in self.gene_directions[field]): #{
      gene_strand = self.gene_directions[field][feature_id]
      if (None == gene_strand): #{
        ExtremeDebugMsg(self, "gene strand not annotated: %s" % feature)
        return
      #} end if
      # if the gene is on the plus strand
      if (gene_strand in "+*"): #{
        strand_flags['plus'] = True
      #} end if
      # if the gene is on the minus strand
      if (gene_strand in "-*"): #{
        strand_flags['minus'] = True
      #} end if
    #} end if
  #} end if

  def OutputEvent(self, event): #{
    if (self.write_special): #{
      mix_count = 0
      for member in event.members: #{
        # check for potential read-through events
        if (member.topology in ["intrachr-same-strand", "read-through"] and
            self.options.read_through > member.meta_fields["target_dist"]): #{
          DebugMsg(self, "Read-through fusion")
          self.outfiles["readthru"].Write(event.FullDataString())
          return
        #} end if
        # if the member has mixed gene directions
        if (hasattr(member, "ends") and None != member.ends and
            (("5","5") == member.ends or ("3","3") == member.ends)): #{
          DebugMsg(self, "Event member with mixed gene directions")
          mix_count += 1
        #} end if
        # if this is an interchromosomal fusion, check whether it
        # involves copies of the same gene on different chromosomes
        if ("interchr" == member.topology and member.GeneSetsOverlap()): #{
          DebugMsg(self, "Interchromosomal fusion involving duplicated gene")
          # write the event to the duplicated-gene output file
          #self.dup_gene_file.Write(event.FullDataString())
          self.outfiles["dup_gene"].Write(event.FullDataString())
          return
        #} end if
      #} end for
      if (len(event.members) == mix_count): #{
        DebugMsg(self, "Event with mixed gene directions")
        # write the event to the mixed-directions output file
        self.outfiles["mix_dirs"].Write(event.FullDataString())
        return
      #} end if
    #} end if
    BasePredictorCls.OutputEvent(self, event)
  #} end def

  #def ReprocessMember(self, member, contig_alignment, contig_seq=None): #{}
  def ReprocessMember(self, member, realigned_contig): #{
    contig_alignment = realigned_contig.best_align
    if (None != contig_alignment and contig_alignment.perfect): #{
      DebugMsg(self, "  False positive: perfect alignment (%i)" %
        contig_alignment.match)
      return False
    #} end if
    # check for poly-A false positives
    align_endA = member.align_info_A.ctg_end
    align_startB = member.align_info_B.ctg_start
    unaligned_end = realigned_contig.sequence[align_endA:].lower()
    unaligned_start = realigned_contig.sequence[:align_startB-1].lower()
    ExtremeDebugMsg(self, "Unaligned End: %s\nUnaligned Start: %s" %
      (unaligned_end, unaligned_start))
    if (unaligned_end.endswith("a"*(len(unaligned_end)-1)) or
        unaligned_start.startswith("t"*(len(unaligned_start)-1))): #{
      DebugMsg(self, "  False positive: poly-A tail (%s or %s)" %
        (unaligned_start, unaligned_end))
      return False
    #} end if
    if (None != contig_alignment): #{
      #new_span = (contig_alignment.ctg_end - contig_alignment.ctg_start) + 1
      #min_span = int(max(member.align_info_A.ContigSpan(),
      #  member.align_info_B.ContigSpan()) * MIN_IMPROVE_FRACTION)
      #DebugMsg(self, "  New span: %i, Min: %i" % (new_span, min_span))
      new_match = contig_alignment.match
      min_match = int(max(member.align_info_A.ContigSpan(),
        member.align_info_B.ContigSpan()) * MIN_IMPROVE_FRACTION)
      DebugMsg(self, "  New match: %i, Min: %i" % (new_match, min_match))
      # if the new alignment is significantly longer (in the contig) than
      # the split contig-to-genome alignments, then the contig is probably
      # not really a fusion
      #if (MIN_IMPROVEMENT < (new_span - member.align_info_A.ContigSpan()) and
      #    MIN_IMPROVEMENT < (new_span - member.align_info_B.ContigSpan())): #{}
      #if (min_span < new_span): #{
      #  DebugMsg(self, "  False positive: better alignment")
      #  return False
      #} end if
      if (min_match < new_match): #{
        DebugMsg(self, "  False positive: better alignment")
        return False
      #} end if
    #} end if
    region_pair = realigned_contig.region_pairs[member.IDString()]
    gene_str_A = "    A (%i-%i): %s" % (region_pair.starts["A"],
      region_pair.ends["A"], ",".join(sorted(region_pair.gene_sets["A"])))
    gene_str_B = "    B (%i-%i): %s" % (region_pair.starts["B"],
      region_pair.ends["B"], ",".join(sorted(region_pair.gene_sets["B"])))
    DebugMsg(self, "\n".join([gene_str_A, gene_str_B]))
    if (region_pair.GeneSetsOverlap()): #{
      DebugMsg(self, "  False positive: gene sets overlap")
      return False
    #} end if
    # check for duplication false positives
    if (self.options.fusion_dup_filter): #{
      DebugMsg(self, "  Applying fusion duplication filter")
      try: #{
        self.missing_align = False
        # check the original contig to genome alignments
        DebugMsg(self, "  Original c2gA")
        if (self.DupFilter(realigned_contig.sequence,
            member.align_info_A, "start")): #{
          DebugMsg(self, "    %s" % member.IDString())
          return False
        #} end if
        DebugMsg(self, "  Original c2gB")
        if (self.DupFilter(realigned_contig.sequence,
            member.align_info_B, "end")): #{
          DebugMsg(self, member.IDString())
          return False
        #} end if
        # check the contig to transcript realignments
        DebugMsg(self, "  c2t realignA")
        if (self.DupFilter(realigned_contig.sequence,
            region_pair.best_aligns['A'], "start")): #{
          DebugMsg(self, member.IDString())
          return False
        #} end if
        DebugMsg(self, "  c2t realignB")
        if (self.DupFilter(realigned_contig.sequence,
            region_pair.best_aligns['B'], "end")): #{
          DebugMsg(self, member.IDString())
          return False
        #} end if
        if (self.missing_align): #{
          self.outfiles["align_miss"].WriteLine("Member %s is missing a "
            "realignment" % member.IDString())
        #} end if
      except AttributeError, e:
        raise FusionPredictorError("Error reprocessing member %s: %s" %
          (member.IDString(), e))
      #} end try
    #} end if
    DebugMsg(self, "  Good fusion prediction")
    return True
  #} end def

  def DupFilter(self, ctg_seq, align, aligned_part): #{
    if (None == align): #{
      self.missing_align = True
      DebugMsg(self, "Missing alignment")
      return False
    #} end if
    ctg_seq = ctg_seq.lower()
    if ("start" == aligned_part): #{
      unaligned_seq = ctg_seq[align.ctg_end:]
    elif ("end" == aligned_part):
      unaligned_seq = ctg_seq[:align.ctg_start-1]
    else:
      raise FusionPredictorError("unrecognized aligned part value: %s" %
        aligned_part)
    #} end if
    aligned_seq = ctg_seq[align.ctg_start-1:align.ctg_end]
    ExtremeDebugMsg(self, "  Unaligned: %s" % unaligned_seq)
    #DebugMsg(self, "Unaligned: %s\nAligned: %s" %
    #  (unaligned_seq, aligned_seq))
    # allow one mismatch at either end
    if ((unaligned_seq[1:] in aligned_seq) or
        (unaligned_seq[:-1] in aligned_seq)): #{
      DebugMsg(self, "  False positive: unaligned in aligned!")
      return True
    #} end if
    if (hasattr(align, "target")): #{
      DebugMsg(self, "    target: %s" % align.target)
      if (align.target in self.transcript_sequences): #{
        DebugMsg(self, "Checking target!")
        transcript_seq = self.transcript_sequences[align.target].lower()
        #DebugMsg(self, "Transcript: %s" % transcript_seq)
        rev_unaligned_seq = ReverseComplement(unaligned_seq)
        ExtremeDebugMsg(self, "  Unaligned_RC: %s\n  Target: %s" %
          (rev_unaligned_seq, align.target))
        if ((unaligned_seq[1:] in transcript_seq) or
            (unaligned_seq[:-1] in transcript_seq) or
            (rev_unaligned_seq[1:] in transcript_seq) or
            (rev_unaligned_seq[:-1] in transcript_seq)): #{
          DebugMsg(self, "  False positive: unaligned in transcript!")
          return True
        #} end if
      else:
        DebugMsg(self, "Target sequence was not preserved:\n%s" %
          ", ".join(sorted(self.transcript_sequences.keys())))
      #} end if
    #} end if
    return False
  #} end def
#} end class

def HandleIntrons(input_list, log_info=None): #{
  ExtremeDebugMsg(log_info, "HANDLE INTRONS")
  if (not isinstance(input_list, MultiDictCls)): #{
    return input_list
  #} end if
  output_list = MultiDictCls("none", input_list.depth,
    list_end=input_list.list_end)
  intron_list = MultiDictCls("none", input_list.depth,
    list_end=input_list.list_end)
  if (2 == input_list.depth): #{
    for gene in input_list: #{
      for part in input_list[gene]: #{
        if (INTRON == part): #{
          intron_list.Add([gene], part)
        else:
          output_list.Add([gene], part)
        #} end if
      #} end for
    #} end for
  elif (4 == input_list.depth):
    for gene in input_list: #{
      for transcript in input_list[gene]: #{
        for part in input_list[gene][transcript]: #{
          keys = [gene, transcript, part]
          for num in input_list[gene][transcript][part]: #{
            if (INTRON == part): #{
              intron_list.Add(keys, num)
            else:
              output_list.Add(keys, num)
            #} end if
          #} end for
        #} end for
      #} end for
    #} end for
  else:
    raise FusionPredictorError("Unrecognize multidictcls depth for "
      "overlapped genes: %i" % input_list.depth)
  #} end if
  ExtremeDebugMsg(log_info, "  NO INTRONS: %s\n  INTRONS: %s" %
    (output_list.OutputString(sort=True, sort_key=str.lower),
     intron_list.OutputString(sort=True, sort_key=str.lower)))
  if (0 == len(output_list)): #{
    return intron_list
  #} end if
  return output_list
#} end def

#### EXCEPTION CLASSES ####
class FusionPredictorError(MyError): #{
  pass
#} end class
