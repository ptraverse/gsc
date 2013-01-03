#! /usr/bin/env python
"""
event_simulator.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, random, re, sys, time, traceback

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
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
  ReverseComplement, NonStandardChr, IntFloor, NormalizeChrID)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  FileBoxCls)
from annotation.create_gene_feature_coords import FixAnnotation
from parsers.genes.annotation import GeneAnnotationParserCls
from parsers.two_bit import TwoBitFileCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"
MAX_COUNTER = 10000

class ChimeraSimulatorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    self.simulator_types = dict()
    self.SimulateChimera = dict()
    if (0 < options.num_fusions): #{
      simulator = SimulatorTypeCls("fusion", options.num_fusions,
        options.output_dir)
      self.simulator_types[simulator.type] = simulator
      self.SimulateChimera[simulator.type] = self.SimulateFusion
    #} end if
    if (0 < options.num_ptds): #{
      simulator = SimulatorTypeCls("ptd", options.num_ptds,
        options.output_dir)
      self.simulator_types[simulator.type] = simulator
      self.SimulateChimera[simulator.type] = self.SimulatePTD
    #} end if
    if (0 < options.num_itds): #{
      simulator = SimulatorTypeCls("itd", options.num_itds,
        options.output_dir)
      self.simulator_types[simulator.type] = simulator
      self.SimulateChimera[simulator.type] = self.SimulateITD
    #} end if
    self.genes_by_chr = dict()
    self.genes_by_id = dict()
    self.counter = 0
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def SimulateChimeras(self): #{
    LogMsg(self, "Simulating chimeric transcripts...")
    start = time.time()
    self.LoadGenes()
    DebugMsg(self, "Loading genome sequence file...")
    genome_file = TwoBitFileCls(self.options.genome_path,
      debug=False) #debug=self.log_info['debug'])
    for simulator in self.simulator_types.itervalues(): #{
      for index in range(simulator.num): #{
        self.ResetCounter()
        chimera = self.SimulateChimera[simulator.type](index,
            genome_file)
        self.GetChimeraSequence(chimera, genome_file)
        if (self.options.min_len < len(chimera)): #{
          simulator.WriteChimera(chimera)
        else:
          ExtremeDebugMsg(self, "Chimera is too short %i (min %i)" %
            (len(chimera), self.options.min_len))
          index -= 1
        #} end if
      #} end for
    #} end for
    LogMsg(self, "Time spent simulating chimeric transcripts: %s" %
      TimeSpent(start))
  #} end def

  def LoadGenes(self): #{
    LogMsg(self, "Loading gene annotations...")
    start = time.time()
    annots_file = AnnotationsFileCls(self.options.annotations_path,
        log_info=self.log_info)
    # iterate through the annotations file
    for transcript in annots_file: #{
      if (None != self.options.chr_filter and
          NormalizeChrID(transcript.chrom) not in self.options.chr_filter): #{
        continue
      #} end if
      # ensure we are not overwriting any other isoforms for
      # the current transcripts
      while (TranscriptID(transcript) in self.genes_by_id): #{
        transcript.isoform += 1
      #} end if
      id = TranscriptID(transcript)
      self.genes_by_id[id] = transcript
      if (transcript.chrom not in self.genes_by_chr): #{
        self.genes_by_chr[transcript.chrom] = dict()
      #} end if
      self.genes_by_chr[transcript.chrom][id] = transcript.alias
    #} end for
    LogMsg(self, "Time spent loading gene annotations: %s" % TimeSpent(start))
  #} end def

  def SimulateFusion(self, fusion_index, genome_file): #{
    DebugMsg(self, "Simulating fusion %i" % fusion_index)
    start = time.time()
    chrB = chrA = self.PickFromList(self.genes_by_chr.keys())
    DebugMsg(self, "Interchromosomal test:")
    if (self.TestValue(self.options.fract_inter)): #{
      DebugMsg(self, "Interchromosomal")
      while (chrB == chrA): #{
        chrB = self.PickFromList(self.genes_by_chr.keys())
      #} end while
    #} end if
    DebugMsg(self, "ChrA: %s, ChrB: %s" % (chrA, chrB))
    gene_idA = self.PickFromList(self.genes_by_chr[chrA].keys())
    while (self.options.skip_single and
           2 > len(self.genes_by_id[gene_idA].exons)): #{
      DebugMsg(self, "Skipping single-exon gene: %s" % gene_idA)
      gene_idA = self.PickFromList(self.genes_by_chr[chrA].keys())
    #} end while
    gene_idB = gene_idA
    while (gene_idB == gene_idA or (self.options.skip_single and
           2 > len(self.genes_by_id[gene_idB].exons))): #{
      gene_idB = self.PickFromList(self.genes_by_chr[chrB].keys())
    #} end while
    geneA = self.genes_by_id[gene_idA]
    geneB = self.genes_by_id[gene_idB]
    DebugMsg(self, "GeneA: %s, GeneB: %s" % (gene_idA, gene_idB))
    # reverse the direction of the first gene half of the time
    DebugMsg(self, "Reverse A test:")
    if (self.TestValue(0.5)): #{
      DebugMsg(self, "Reversing direction of Gene A: seq will be 3'<5'")
      reverseA = True
    else:
      DebugMsg(self, "Not reversing direction of Gene A: seq will be 5'>3'")
      reverseA = False
    #} end if
    # determine whether transcription direction should be maintained
    DebugMsg(self, "Maintain direction test:")
    if (self.TestValue(self.options.fract_cont)):
      if (reverseA): #{
        DebugMsg(self, "Maintaining direction by reversing B: "
          "seq will be 3'<5'")
        reverseB = True
      else:
        DebugMsg(self, "Maintaining direction by not reversing B: "
          "seq will be 5'>3'")
        reverseB = False
      #} end if
    else:
      if (reverseA): #{
        DebugMsg(self, "Not maintaining direction by not reversing B: "
          "seq will be 5'>3'")
        reverseB = False
      else:
        DebugMsg(self, "Not maintaining direction by reversing B: "
          "seq will be 3'<5'")
        reverseB = True
      #} end if
    #} end if
    # determine whether to restrict breakpoints to exon-edges
    DebugMsg(self, "Exon breakpoint test:")
    if (self.TestValue(self.options.fract_exon) and
        1 < len(geneA.exons) and 1 < len(geneB.exons)): #{
      DebugMsg(self, "Randomly choosing exons as breakpoints")
      use_exons = True
      splitA = self.PickValue(len(geneA.exons)-2)
      splitB = self.PickValue(len(geneB.exons)-2)
    else:
      if (1 == len(geneA.exons) or 1 == len(geneB.exons)): #{
        DebugMsg(self, "Single-exon gene, defaulting to position")
      #} end if
      DebugMsg(self, "Randomly choosing positions as breakpoints")
      use_exons = False
      splitA = self.RandomTranscriptPosition(geneA)
      splitB = self.RandomTranscriptPosition(geneB)
    #} end if
    DebugMsg(self, "SplitA: %i, SplitB: %i" % (splitA, splitB))
    DebugMsg(self, "Creating fusion transcript from:\n"
      "  GeneA: %s; %s\n" % (geneA.details(), geneA.exon_string()) +
      "  GeneB: %s; %s" % (geneB.details(), geneB.exon_string()))
    fusion = FusionTranscriptCls(fusion_index, geneA, geneB,
      reverseA, reverseB, use_exons, splitA, splitB)
    DebugMsg(self, "Fusion: %s\n" % fusion.ToString() +
      "Time spent simulating fusion: %s" % TimeSpent(start))
    return fusion
  #} end def

  def SimulatePTD(self, ptd_index, genome_file): #{
    DebugMsg(self, "Simulating PTD %i" % ptd_index)
    start = time.time()
    gene_id = self.PickFromList(self.genes_by_id.keys())
    gene = self.genes_by_id[gene_id]
    while (3 > len(gene.exons) or
           (3 == len(gene.exons) and 2 > ExonLength(gene.exons[1]))): #{
      if (3 == len(gene.exons)): #{
        DebugMsg(self, "Skipping gene with only single-base internal exon")
      else:
        DebugMsg(self, "Skipping gene with too few exons: %s" % gene_id)
      #} end if
      gene_id = self.PickFromList(self.genes_by_id.keys())
      gene = self.genes_by_id[gene_id]
    #} end while
    DebugMsg(self, "Gene: %s" % (gene_id))
    exonA = self.PickInternalExon(gene)
    while (2 > ExonLength(gene.exons[exonA])): #{
      if (2 > ExonLength(gene.exons[exonA])): #{
        DebugMsg(self, "Skipping single-base exon: %i" % exonA)
      else:
        DebugMsg(self, "Skipping UTR: %i" % exonA)
      #} end if
      exonA = self.PickInternalExon(gene)
    #} end while
    DebugMsg(self, "Multiple exons test")
    exonB = exonA
    if (3 < len(gene.exons) and
        self.TestValue(self.options.fract_multi)): #{
      DebugMsg(self, "Selecting second exon")
      while (exonB == exonA): #{
        exonB = self.PickInternalExon(gene)
      #} end while
    #} end if
    DebugMsg(self, "ExonA: %i, ExonB: %i" % (exonA, exonB))
    extra_seq = ""
    DebugMsg(self, "Extra sequence test")
    if (self.TestValue(self.options.fract_extra_ptd)): #{
      DebugMsg(self, "Generating extra sequence")
      extra_seq = self.RandomSequence(self.options.min_extra_ptd,
        self.options.max_extra_ptd)
    #} end if
    DebugMsg(self, "Creating PTD transcript from:\n"
      "  Gene: %s; %s" % (gene.details(), gene.exon_string()))
    ptd = PTDTranscriptCls(ptd_index, gene, exonA, exonB, extra_seq)
    DebugMsg(self, "PTD: %s\n" % ptd.ToString() +
      "Time spent simulating PTD: %s" % TimeSpent(start))
    return ptd
  #} end def

  def SimulateITD(self, itd_index, genome_file): #{
    DebugMsg(self, "Simulating ITD %i" % itd_index)
    start = time.time()
    gene_id = self.PickFromList(self.genes_by_id.keys())
    gene = self.genes_by_id[gene_id]
    max_exon_len = GetMaxExonLength(gene)
    while ((self.options.skip_single and 2 > gene.num_coding_exons) or
           (1.5*self.options.min_len_itd > max_exon_len) or
           (gene.non_coding)): #{
      if (1.5*self.options.min_len_itd > max_exon_len): #{
        DebugMsg(self, "Skipping gene with all short exons: %s" % gene_id)
      elif (gene.non_coding):
        DebugMsg(self, "Skipping non-coding gene: %s" % gene_id)
      else:
        DebugMsg(self, "Skipping single-exon gene: %s" % gene_id)
      #} end if
      gene_id = self.PickFromList(self.genes_by_id.keys())
      gene = self.genes_by_id[gene_id]
      max_exon_len = GetMaxExonLength(gene)
    #} end while
    DebugMsg(self, "Gene: %s Max exon: %i" % (gene_id, max_exon_len))
    if ("+" == gene.strand): #{
      include_pos = False
    else:
      include_pos = True
    #} end if
    # pick an exon (ensure it is large enough)
    exon = self.PickValue(len(gene.split_exons)-1)
    exon_length = ExonLength(gene.split_exons[exon])
    while (1.5*self.options.min_len_itd > exon_length or
           gene.utr_flags[exon]): #{
      if (gene.utr_flags[exon]): #{
        DebugMsg(self, "Skipping UTR: %i" % exon)
      else:
        DebugMsg(self, "Skipping too-short exon: %i (%i)" %
            (exon, exon_length))
      #} end if
      exon = self.PickValue(len(gene.split_exons)-1)
      exon_length = ExonLength(gene.split_exons[exon])
    #} end while
    DebugMsg(self, "Exon: %s, Length: %s" % (exon, exon_length))
    # get the start and end positions (ensure size is valid)
    posA = self.PickValue(exon_length-1)
    posB = posA
    dist = abs(posA-posB)+1
    dup_seq = self.GetDupSeq(gene, exon, posA, posB, genome_file)
    while (self.options.min_len_itd > dist or
           self.options.max_len_itd < dist or
           self.HomopolymerTest(dup_seq)): #{
      posB = self.PickValue(exon_length-1)
      dist = abs(posA-posB)+1
      dup_seq = self.GetDupSeq(gene, exon, posA, posB, genome_file)
      DebugMsg(self, "Duplication size: %i Valid: %i-%i" %
          (dist, self.options.min_len_itd, self.options.max_len_itd))
    #} end while
    DebugMsg(self, "PosA: %i, PosB: %i" % (posA, posB))
    # get extra sequence
    extra_seq = ""
    DebugMsg(self, "Extra sequence test")
    if (self.TestValue(self.options.fract_extra_itd)): #{
      DebugMsg(self, "Generating extra sequence")
      extra_seq = self.RandomSequence(self.options.min_extra_itd,
        self.options.max_extra_itd)
    #} end if
    DebugMsg(self, "Creating ITD transcript from:\n"
      "  Gene: %s; %s" % (gene.details(), gene.exon_string()))
    itd = ITDTranscriptCls(itd_index, gene, exon, posA, posB, extra_seq)
    DebugMsg(self, "ITD: %s\n" % itd.ToString() +
      "Time spent simulating ITD: %s" % TimeSpent(start))
    return itd
  #} end def

  def GetDupSeq(self, gene, exon, posA, posB, genome_file): #{
    exon = gene.split_exons[exon]
    dup_start = exon[0] + min(posA, posB)
    dup_end   = exon[0] + max(posA, posB)
    return genome_file.GetSequence(gene.chrom, dup_start, dup_end)
  #} end def

  def HomopolymerTest(self, test_seq): #{
    homopolymer_patt = r"^(.)\1*$"
    if (None == re.search(homopolymer_patt, test_seq, re.IGNORECASE)): #{
      return False
    #} end if
    return True
  #} end def

  def ResetCounter(self): #{
    self.counter = 0
  #} end def

  def UpdateCounter(self): #{
    self.counter += 1
    if (self.counter > MAX_COUNTER): #{
      raise ChimeraSimulatorError("loop overflow while simulating "
        "chimeric transcript")
    #} end if
  #} end def

  def TestValue(self, test_value): #{
    rand_value = random.random()
    DebugMsg(self, "Random: %f, Test: %f" % (rand_value, test_value))
    return rand_value < test_value
  #} end def

  def PickInternalExon(self, gene): #{
    return self.PickValue(len(gene.exons)-3) + 1
  #} end def

  def PickFromList(self, source_list): #{
    self.UpdateCounter()
    return random.choice(source_list)
  #} end def

  def PickValue(self, max_value): #{
    self.UpdateCounter()
    return random.randint(0, max_value)
  #} end def

  def RandomTranscriptPosition(self, transcript): #{
    length = TranscriptLength(transcript)
    DebugMsg(self, "Length: %i" % (length))
    max_value = length - 2*self.options.min_len_fus
    return self.options.min_len_fus + self.PickValue(max_value)
  #} end def

  def RandomSequence(self, min_len, max_len): #{
    bases = "ACGT"
    sequence = ""
    length = random.randint(min_len, max_len)
    while (len(sequence) < length): #{
      sequence += random.choice(bases)
    #} end while
    DebugMsg(self, "Length: %i, Seq: %s" % (length, sequence))
    return sequence
  #} end def

  def GetChimeraSequence(self, chimera, genome_file): #{
    # clear the current sequence
    chimera.seq = ""
    chimera.len = 0
    # add sequence for each block in the first part of the chimera
    DebugMsg(self, "Getting sequence for blocks before event")
    self.GetSequencePart(chimera, genome_file, chimera.blocksA)
    if (None != chimera.extra_seq and "" != chimera.extra_seq): #{
      DebugMsg(self, "Adding extra sequence: %i" % len(chimera.extra_seq))
      chimera.seq += chimera.extra_seq
      chimera.len += len(chimera.extra_seq)
    #} end if
    # add sequence for each block in the second part of the chimera
    DebugMsg(self, "Getting sequence for blocks after event")
    self.GetSequencePart(chimera, genome_file, chimera.blocksB)
    DebugMsg(self, "Chimera size: %i, Sequence length: %i" %
      (chimera.len, len(chimera.seq)))
    if (chimera.len != len(chimera.seq)): #{
      raise ChimeraSimulatorError("incorrect sequence length!")
    #} end if
  #} end def

  def GetSequencePart(self, chimera, genome_file, block_list): #{
    for block in block_list: #{
      new_seq = genome_file.GetSequence(block.chr, block.start, block.end)
      if (block.apply_rc): #{
        new_seq = ReverseComplement(new_seq)
      #} end if
      chimera.seq += new_seq
      chimera.len += block.Length()
      #DebugMsg(self, "Block size: %i, Sequence Length: %i" %
      #  (block.Length(), len(new_seq)))
    #} end for
  #} end def
#} end class

class SimulatorTypeCls: #{
  def __init__(self, type, num, output_dir): #{
    self.type = type
    self.num = num
    output_path = os.path.join(output_dir, "simulated_%ss.fa" % type)
    self.file = FileBoxCls(output_path, "w",
      "cannot create simulated %s output path")
  #} end def

  def WriteChimera(self, chimera): #{
    self.file.WriteLine(">%s" % chimera.ToString())
    self.file.WriteLine(chimera.seq)
  #} end def
#} end class

class AnnotationsFileCls: #{
  def __init__(self, path, log_info=None): #{
    self.log_info = log_info
    # open the annotations file
    self.parser = GeneAnnotationParserCls(path, log_info=log_info)
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  def next(self): #{
    # replace spaces in the alias and
    # get rid of any "chr" in the chromosome name
    transcript = FixAnnotation(self.parser.next(), use_chr=False)
    ExtremeDebugMsg(self, "T: %s" % transcript)
    # ensure that the transcript is from a "normal" chromosome,
    # not including mitochondrial DNA, and is not a tRNA or rRNA
    while (NonStandardChr(transcript.chrom) or
        "M" == transcript.chrom or
        transcript.gene_name.lower().startswith("trna_") or
        transcript.gene_name.lower().endswith("_rrna")): #{
      ExtremeDebugMsg(self, "  Skipping...")
      transcript = FixAnnotation(self.parser.next(), use_chr=False)
      ExtremeDebugMsg(self, "T: %s" % transcript)
    #} end while
    transcript.isoform = 1
    # check whether the transcript is coding or non-coding
    if (transcript.cdsStart >= transcript.cdsEnd): #{
      transcript.non_coding = True
    #} end if
    # separate the exons into UTRs and coding exons
    self.SeparateUTRs(transcript)
    # reverse the order of the exons if
    # the transcript is on the negative strand
    if ("-" == transcript.strand): #{
      transcript.exons.reverse()
      transcript.split_exons.reverse()
      transcript.utr_flags.reverse()
    #} end if
    return transcript
  #} end def

  def SeparateUTRs(self, transcript): #{
    # ensure that exons are ordered by start coordinate
    if (transcript.exons[0][0] > transcript.exons[-1][0]): #{
      transcript.exons.reverse()
    #} end if
    transcript.num_coding_exons = 0
    transcript.utr_flags = list()
    transcript.split_exons = list()
    if (transcript.non_coding): #{
      ExtremeDebugMsg(self, "Not separating UTRs for non-coding gene")
      return
    #} end if
    ExtremeDebugMsg(self, "Separating UTRs from coding exons...\n"
      "cdsStart: %i, cdsEnd: %i" % (transcript.cdsStart, transcript.cdsEnd))
    for (e_start, e_end) in transcript.exons: #{
      ExtremeDebugMsg(self, "Exon start: %i, end: %i" % (e_start, e_end))
      # if the exon ends before the CDS start or
      # the exon starts after the CDS end,
      # the full exon is a UTR
      if (e_end < transcript.cdsStart or transcript.cdsEnd < e_start): #{
        transcript.utr_flags.append(True)
        transcript.split_exons.append([e_start, e_end])
        ExtremeDebugMsg(self, "  full UTR")
      else:
        # if the exon starts before the CDS start and
        # ends after the CDS start,
        # the first part of the exon is a UTR
        if (e_start < transcript.cdsStart): #{
          transcript.utr_flags.append(True)
          transcript.split_exons.append([e_start, transcript.cdsStart-1])
          e_start = transcript.cdsStart
          ExtremeDebugMsg(self, "  UTR start: %i-%i\n  New start: %i" %
            (transcript.split_exons[-1][0],
            transcript.split_exons[-1][1], e_start))
        #} end if
        # if the exon starts before the CDS end and
        # ends after the CDS end,
        # the second part of the exon is a UTR
        if (transcript.cdsEnd < e_end): #{
          transcript.num_coding_exons += 1
          transcript.utr_flags.append(False)
          transcript.split_exons.append([e_start, transcript.cdsEnd])
          transcript.utr_flags.append(True)
          transcript.split_exons.append([transcript.cdsEnd+1, e_end])
          ExtremeDebugMsg(self, "  exon start: %i-%i\n  UTR end: %i-%i" %
            (transcript.split_exons[-2][0], transcript.split_exons[-2][1],
             transcript.split_exons[-1][0], transcript.split_exons[-1][1]))
        # if the exon starts after the CDS start and
        # ends before the CDS end,
        # the full exon is really an exon
        elif (e_start <= e_end):
          transcript.num_coding_exons += 1
          transcript.utr_flags.append(False)
          transcript.split_exons.append([e_start, e_end])
          ExtremeDebugMsg(self, "  full exon: %i-%i" % (e_start, e_end))
        else:
          raise ExonCoordsError("cannot determine exon type: "
            "%s: CDS:%i-%i, Exon:%i-%i" % (transcript.alias,
            transcript.cdsStart, transcript.cdsEnd, e_start, e_end))
        #} end if
      #} end if
    #} end for
    if (len(transcript.split_exons) != len(transcript.utr_flags)): #{
      raise ChimeraSimulatorError("error loading transcript: # exons (%i)" %
        len(transcript.exons) + " not equal to # UTR flags (%i)" %
        len(transcript.utr_flags))
    #} end if
  #} end def
#} end class

class ChimeraCls: #{
  def __init__(self, type, index): #{
    self.type      = type
    self.index     = index
    self.gene      = None
    self.blocksA   = list()
    self.blocksB   = list()
    self.extra_seq = None
    self.seq       = ""
  #} end def

  def __len__(self): #{
    return len(self.seq)
  #} end def

  def ToString(self): #{
    split = sum(map(len, self.blocksA))
    data_list = [
      "%s%i" % (self.type, self.index),
      "SPLIT:%i" % split,
      "GENE:%s" % self.gene,
      self.TypeSpecificString(),
    ]
    if (None != self.extra_seq and "" != self.extra_seq): #{
      data_list.append("EXTRA:%s" % self.extra_seq)
    #} end if
    block_strA = ",".join(block.ToString() for block in self.blocksA)
    block_strB = ",".join(block.ToString() for block in self.blocksB)
    data_list.append("BLOCKS:%s" % "/".join([block_strA, block_strB]))
    return " ".join(data_list)
  #} end def
#} end class

class FusionTranscriptCls(ChimeraCls): #{
  def __init__(self, index, geneA, geneB, reverseA, reverseB,
    use_exons, splitA, splitB): #{
    ChimeraCls.__init__(self, 'fusion', index)
    self.gene = "%s(%s)/%s(%s)" % (geneA.alias, geneA.name,
      geneB.alias, geneB.name)
    if (reverseA): #{
      strandA = OppositeStrand(geneA.strand)
    else:
      strandA = geneA.strand
    #} end if
    if (reverseB): #{
      strandB = OppositeStrand(geneB.strand)
    else:
      strandB = geneB.strand
    #} end if
    self.chr = "chr%s(%s)/chr%s(%s)" % (geneA.chrom, strandA,
      geneB.chrom, strandB)
    self.blocksA = GetBlocksFromGene(geneA, reverseA, use_exons, splitA)[0]
    self.blocksB = GetBlocksFromGene(geneB, reverseB, use_exons, splitB)[1]
    self.dirA = GetDirectionString(reverseA)
    self.dirB = GetDirectionString(reverseB)
    self.use_exons = use_exons
    #if (reverseA): #{
    #  include_posA = True
    #else:
    #  include_posA = False
    #} end if
    if (reverseA): #{
      self.splitA = splitA+1
    else:
      self.splitA = splitA
    #} end if
    #self.splitA = GetDirectedValue(geneA, splitA, use_exons, include_posA)
    #if (reverseB): #{
    #  include_posB = False
    #else:
    #  include_posB = True
    #} end if
    if (reverseB): #{
      self.splitB = splitB
    else:
      self.splitB = splitB+1
    # end if
    #self.splitB = GetDirectedValue(geneB, splitB, use_exons, include_posB)
  #} end def

  def TypeSpecificString(self): #{
    data_list = [
      "CHRS:%s" % self.chr,
      "DIRECTIONS:%s/%s" % (self.dirA, self.dirB),
      "BREAKPOINTS:%s/%s" %
        (self.BreakPoint(self.splitA), self.BreakPoint(self.splitB)),
    ]
    return " ".join(data_list)
  #} end def

  def BreakPoint(self, index): #{
    if (self.use_exons): #{
      return "exon%i" % (index+1)
    else:
      return "pos%s" % index
    #} end if
  #} end def
#} end class

class PTDTranscriptCls(ChimeraCls): #{
  def __init__(self, index, gene, exonA, exonB, extra_seq): #{
    ChimeraCls.__init__(self, 'ptd', index)
    self.gene = "%s(%s)" % (gene.alias, gene.name)
    self.dup_start = min(exonA, exonB)
    self.dup_end   = max(exonA, exonB)
    self.dup_len   = 0
    for exon in gene.exons[self.dup_start:self.dup_end+1]: #{
      self.dup_len += ExonLength(exon)
    #} end for
    self.total_len = self.dup_len + len(extra_seq)
    #if ("+" == gene.strand): #{
    #  include_pos = False
    #  dup_start = min(exonA, exonB)-1
    #  dup_end   = max(exonA, exonB)
    #else:
    #  include_pos = True
    #  dup_start = max(exonA, exonB)
    #  dup_end   = min(exonA, exonB)-1
    #} end if
    self.blocksA = GetBlocksFromGene(gene, False, True, self.dup_end)[0]
    self.blocksB = GetBlocksFromGene(gene, False, True, self.dup_start-1)[1]
    #directedA = GetDirectedValue(gene, exonA, True, include_pos)
    #directedB = GetDirectedValue(gene, exonB, True, include_pos)
    #self.dup_start = min(directedA, directedB)
    #self.dup_end   = max(directedA, directedB)
    self.extra_seq = extra_seq
  #} end def

  def TypeSpecificString(self): #{
    if (self.dup_start != self.dup_end): #{
      dup_string = "DUP:exon%i-exon%i" % (self.dup_start+1, self.dup_end+1)
    else:
      dup_string = "DUP:exon%i" % (self.dup_start+1)
    #} end if
    return "%s DUP_LEN:%i(%i) TOTAL_LEN:%i(%i)" % (dup_string, self.dup_len,
      (self.dup_len % 3), self.total_len, (self.total_len % 3))
  #} end def
#} end class

class ITDTranscriptCls(ChimeraCls): #{
  def __init__(self, index, gene, exon, posA, posB, extra_seq): #{
    ChimeraCls.__init__(self, 'itd', index)
    self.gene = "%s(%s)" % (gene.alias, gene.name)
    splitA = exon-1
    splitB = exon
    #if ("+" == gene.strand): #{
    #  splitA = exon-1
    #  splitB = exon
    #else:
    #  splitA = exon
    #  splitB = exon-1
    #} end if
    self.blocksA = GetBlocksFromGene(gene, False, True, splitA,
        split_exons=True)[0]
    self.blocksB = GetBlocksFromGene(gene, False, True, splitB,
        split_exons=True)[1]
    #if ("+" == gene.strand): #{
    #  include_pos = False
    #else:
    #  include_pos = True
    #} end if
    self.exon_id = exon
    #self.exon_id = GetDirectedValue(gene, exon, True, include_pos)
    self.SplitExonWithDup(gene, posA, posB)
    self.dup_len = abs(posA - posB) + 1
    self.total_len = self.dup_len + len(extra_seq)
    self.extra_seq = extra_seq
  #} end def

  def SplitExonWithDup(self, gene, posA, posB): #{
    exon = gene.split_exons[self.exon_id]
    #ErrMsg("DEBUG: exon length: %i" % ExonLength(exon))
    start = min(posA, posB)
    end   = max(posA, posB)
    if ("+" == gene.strand): #{
      self.dup_start = start
      self.dup_end   = end
      before_gap = BlockCls(gene.chrom, exon[0], exon[0]+end, False)
      after_gap  = BlockCls(gene.chrom, exon[0]+start, exon[1], False)
    else:
      self.dup_start = ExonLength(exon) - end   - 1
      self.dup_end   = ExonLength(exon) - start - 1
      before_gap = BlockCls(gene.chrom, exon[0]+start, exon[1], True)
      after_gap  = BlockCls(gene.chrom, exon[0], exon[0]+end, True)
    #} end if
    self.blocksA.append(before_gap)
    self.blocksB.insert(0, after_gap)
  #} end def

  def TypeSpecificString(self): #{
    return "EXON:%i DUP:%i-%i DUP_LEN:%i(%i) TOTAL_LEN:%i(%i)" % (
      self.exon_id+1, self.dup_start+1, self.dup_end+1,
      self.dup_len, (self.dup_len % 3), self.total_len, (self.total_len % 3))
  #} end def
#} end class

class BlockCls: #{
  def __init__(self, chr, start, end, apply_rc): #{
    self.chr      = chr
    self.start    = start
    self.end      = end
    self.apply_rc = apply_rc
  #} end def

  def __len__(self): #{
    return self.Length()
  #} end def

  def Length(self): #{
    return (self.end - self.start) + 1
  #} end def

  def ToString(self): #{
    if (self.apply_rc): #{
      return "%sr:%i-%i" % (self.chr, self.end, self.start)
    else:
      return "%s:%i-%i" % (self.chr, self.start, self.end)
    #} end if
  #} end def
#} end class

def TranscriptID(transcript): #{
  return "%s.i%i" % (transcript.alias.lower(), transcript.isoform)
#} end def

def TranscriptLength(transcript): #{
  if (not hasattr(transcript, "length") or 0 == transcript.length): #{
    transcript.length = 0
    for exon in transcript.exons: #{
      transcript.length += ExonLength(exon)
    #} end for
  #} end if
  return transcript.length
#} end def

def ExonLength(exon): #{
  return (exon[1] - exon[0]) + 1
#} end def

def GetMaxExonLength(gene): #{
  max_length = 0
  for exon in gene.split_exons: #{
    exon_length = ExonLength(exon)
    max_length = max(exon_length, max_length)
  #} end for
  return max_length
#} end def

def GetBlocksFromGene(gene, reverse, use_exons, split, split_exons=False): #{
  apply_rc = False
  if (("+" == gene.strand and reverse) or
      ("-" == gene.strand and not reverse)): #{
    apply_rc = True
    #if (use_exons): #{
    #  split = len(gene.exons) - split - 2
    #else:
    #if (not use_exons): #{
    #  split = TranscriptLength(gene) - split - 2
    #} end if
  #} end if
  if (split_exons): #{
    exon_list = gene.split_exons
  else:
    exon_list = gene.exons
  #} end if
  before_split = list()
  after_split = list()
  if (use_exons): #{
    #ErrMsg("Split: %i" % split)
    for (index, exon) in enumerate(exon_list): #{
      #ErrMsg("Processing block %i" % index)
      new_block = BlockCls(gene.chrom, exon[0], exon[1], apply_rc)
      if (index <= split): #{
        #ErrMsg("  Before split")
        before_split.append(new_block)
      else:
        #ErrMsg("  After split")
        after_split.append(new_block)
      #} end if
    #} end for
  else:
    used = 0
    for exon in exon_list: #{
      new_block = BlockCls(gene.chrom, exon[0], exon[1], apply_rc)
      if ((used + new_block.Length()) <= split): #{
        before_split.append(new_block)
      elif (used < split): #{
        remaining = split - used
        if (apply_rc): #{
          before_block = BlockCls(gene.chrom, (exon[1]-remaining)+1, exon[1],
            apply_rc)
          after_block = BlockCls(gene.chrom, exon[0], (exon[1]-remaining),
            apply_rc)
        else:
          before_block = BlockCls(gene.chrom, exon[0], (exon[0]+remaining)-1,
            apply_rc)
          after_block = BlockCls(gene.chrom, (exon[0]+remaining), exon[1],
            apply_rc)
        #} end if
        before_split.append(before_block)
        after_split.append(after_block)
      else:
        after_split.append(new_block)
      #} end if
      used += new_block.Length()
    #} end for
  #} end if
  if (reverse): #{
    before_split.reverse()
    after_split.reverse()
    return (after_split, before_split)
  #} end if
  return (before_split, after_split)
#} end def

def OppositeStrand(in_strand): #{
  if ("+" == in_strand): #{
    return "-"
  elif ("-" == in_strand):
    return "+"
  #} end if
  raise ChimeraSimulatorError("unrecognized strand: %s" % in_strand)
#} end def

def GetDirectionString(reverse): #{
  if (reverse): #{
    return "3'<5'"
  #} end if
  return "5'>3'"
#} end def

def GetDirectedValue(gene, value, use_exons, include_pos): #{
  if ("+" == gene.strand): #{
    if (include_pos): #{
      return value+1
    #} end if
    return value
  #} end if
  if (use_exons): #{
    value = len(gene.exons) - value - 1
  else:
    value = TranscriptLength(gene) - value - 1
  #} end if
  if (include_pos): #{
    return value
  #} end if
  return value - 1
#} end def

#### EXCEPTION CLASSES ####
class ChimeraSimulatorError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = "Produce sequences for simulated chimeric transcripts"
  args = [ "GENE_ANNOTATIONS_FILE", "GENOME_SEQUENCE_FILE", "OUTPUT_DIR", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--skip-single-exon-genes",
                    action="store_true", dest="skip_single",
                    help="Do not use genes with only a single exon when "
                         "simulating events. [default]")
  parser.add_option("--allow-single-exon-genes",
                    action="store_false", dest="skip_single",
                    help="Allow the use of genes with only a single exon "
                         "when simulating events.")
  parser.add_option("--chr-filter",
      help="Only simulate events involving chromosomes in the provided "
        "comma-delimited list.")
  parser.add_option("--min-len",
      type="int", metavar="N",
      help="Ensure that each simulated transcript is at least Nbp long. "
        "[default: %default]")
  fusion_grp = OptionGroup(parser, "Fusion Options")
  fusion_grp.add_option("--num-fusions",
                        type="int", metavar="N",
                        help="Simulate N fusion transcripts. "
                             "[default: %default]")
  fusion_grp.add_option("--fract-inter",
                        type="float", metavar="F",
                        help="The fraction of fusions to be interchromosomal "
                             "(between genes on distinct chromosomes). "
                             "[default: %default]")
  fusion_grp.add_option("--fract-cont",
                        type="float", metavar="F",
                        help="The fraction of fusions to maintain continuous "
                             "transcription direction of both genes. "
                             "[default: %default]")
  fusion_grp.add_option("--fract-exon",
                        type="float", metavar="F",
                        help="The fraction of fusions to have their "
                             "breakpoints exactly at exon boundaries. "
                             "[default: %default]")
  fusion_grp.add_option("--min-len-fus",
                        type="int", metavar="N",
                        help="When the breakpoint is not being constrained "
                             "to exon edges, it must be at least Nbp away "
                             "from the edge of the transcript. "
                             "[default: %default]")
  parser.add_option_group(fusion_grp)
  ptd_grp = OptionGroup(parser, "Partial Tandem Duplication Options")
  ptd_grp.add_option("--num-ptds",
                     type="int", metavar="N",
                     help="Simulate N partial tandem duplication transcripts. "
                          "[default: %default]")
  ptd_grp.add_option("--fract-multi",
                     type="float", metavar="F",
                     help="The fraction of PTDs to involve a duplication "
                          "of multiple exons, rather than just one. "
                          "[default: %default]")
  ptd_grp.add_option("--fract-extra-ptd",
                     type="float", metavar="F",
                     help="The fraction of PTDs to involve extra sequence "
                          "between the copies of the duplicated sequence. "
                          "[default: %default]")
  ptd_grp.add_option("--min-extra-ptd",
                     type="int", metavar="N",
                     help="When extra sequence is inserted in PTDs, it is at "
                          "least Nbp long. [default: %default]")
  ptd_grp.add_option("--max-extra-ptd",
                     type="int", metavar="N",
                     help="When extra sequence is inserted in PTDs, it is at "
                          "most Nbp long. [default: %default]")
  parser.add_option_group(ptd_grp)
  itd_grp = OptionGroup(parser, "Internal Tandem Duplication Options")
  itd_grp.add_option("--num-itds",
                     type="int", metavar="N",
                     help="Simulate N internal tandem duplication "
                          "transcripts. [default: %default]")
  itd_grp.add_option("--min-len-itd",
                     type="int", metavar="N",
                     help="For ITDs, the duplicated sequence must be at least "
                          "Nbp long. [default: %default]")
  itd_grp.add_option("--max-len-itd",
                     type="int", metavar="N",
                     help="For ITDs, the duplicated sequence must be at most "
                          "Nbp long. [default: %default]")
  itd_grp.add_option("--fract-extra-itd",
                     type="float", metavar="F",
                     help="The fraction of ITDs to involve extra sequence "
                          "between the copies of the duplicated sequence. "
                          "[default: %default]")
  itd_grp.add_option("--min-extra-itd",
                     type="int", metavar="N",
                     help="When extra sequence is inserted in ITDs, it is at "
                          "least Nbp long. [default: %default]")
  itd_grp.add_option("--max-extra-itd",
                     type="int", metavar="N",
                     help="When extra sequence is inserted in ITDs, it is at "
                          "most Nbp long. [default: %default]")
  parser.add_option_group(itd_grp)
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--seed",
      type="float",
      help="The seed to use to initialize the random number generator. If no "
        "value is specified, current system time is used.")
  misc_group.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(skip_single=True,
                      min_len=200,
                      num_fusions=0,
                      fract_inter=.5,
                      fract_cont=.75,
                      fract_exon=.9,
                      min_len_fus=15,
                      num_ptds=0,
                      fract_multi=.75,
                      fract_extra_ptd=.05,
                      min_extra_ptd=1,
                      max_extra_ptd=5,
                      num_itds=0,
                      min_len_itd=5,
                      max_len_itd=100,
                      fract_extra_itd=.5,
                      min_extra_itd=1,
                      max_extra_itd=20,
                      force=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (0 > options.num_fusions): #{
    ErrMsg("Invalid number of fusions: %i" % options.num_fusions)
    opts_good = False
  #} end if
  if (0 > options.num_ptds): #{
    ErrMsg("Invalid number of partial tandem duplications: %i" %
      options.num_ptds)
    opts_good = False
  #} end if
  if (0 > options.num_itds): #{
    ErrMsg("Invalid number of internal tandem duplications: %i" %
      options.num_itds)
    opts_good = False
  #} end if
  if (0 > (options.num_fusions + options.num_ptds, options.num_itds)): #{
    ErrMsg("You have not asked to simulate any events!")
    opts_good = False
  #} end if
  path_errors = list()
  # check the gene annotations file path
  CheckFilePath(options.annotations_path, "gene annotation", path_errors)
  # check the genome sequence file path
  CheckFilePath(options.genome_path, "genome sequence", path_errors)
  # check the output path
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.annotations_path,
      "simulator", options.output_dir)
  #} end if
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # seed the random number generator
  if (None == options.seed): #{
    options.seed = IntFloor(time.time())
  #} end if
  random.seed(options.seed)
  # set up the chromosome filter list, if one has been given
  if (None != options.chr_filter): #{
    options.chr_filter = set(map(NormalizeChrID,
      options.chr_filter.split(",")))
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
    options.annotations_path = EnsureAbsPath(args[0])
    options.genome_path      = EnsureAbsPath(args[1])
    options.output_dir       = EnsureAbsPath(args[2])
    if (CheckPaths(options)): #{
      try: #{
        chimera_simulator = ChimeraSimulatorCls(options)
        WriteCommand(chimera_simulator, sys.argv)
        chimera_simulator.SimulateChimeras()
      except (MyError), e:
        ErrMsg("ERROR while simulating chimeric transcripts:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify the path to a gene annotations file "
      "(GENE_ANNOTATIONS_FILE); the path to a 2-bit per base formatted genome "
      "sequence file (GENOME_SEQUENCE_FILE); and the path to an output "
      "directory (OUTPUT_DIR).")
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
