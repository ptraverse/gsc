"""
transcript.py

Created by Readman Chiu
Edited by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

import os, re
from parsers import genes
from utils.general import IsInteger
from common.coord_pair import CoordPairCls

class Transcript:

  # name = transcript ID
  name = None
  chrom = None
  strand = None
  txStart = txEnd = None
  cdsStart = cdsEnd = None
  exonCount = None
  exons = None
  gene = None
  # alias = gene name
  alias = None
  model = None
  weight = None
  length = 0

  def __init__(self, name):
    self.name = name

  def full_name(self):
    alias = "NA"
    if self.alias:
      alias = self.alias
    return "%s(%s)" % (self.name, alias)

  def details(self):
    return "%s %s %s %s-%s %s" % (self.full_name(), self.chrom, self.strand, self.txStart, self.txEnd, self.exonCount)

  def exon_string(self):
    return ",".join(["%i-%i" % (exon[0], exon[1]) for exon in self.exons])

  def cds_length(self):
    if self.cdsStart and self.cdsEnd and IsInteger(self.cdsStart) and IsInteger(self.cdsEnd):
      return int(self.cdsEnd) - int(self.cdsStart)
    return 0

  def coding_type(self):
    if self.cdsStart == None or self.cdsEnd == None:
      return 'NA'
    elif IsInteger(self.cdsStart) and IsInteger(self.cdsEnd) and self.cdsStart != self.cdsEnd:
      return 'CODING'
    else:
      return 'NONCODING'

  @classmethod
  def find_by_name(cls, model, annot_file, txt_name):
    txts = []
    if os.path.exists(annot_file):
      for line in open(annot_file, 'r'):
        if txt_name in line:
          txt = {
            'e': genes.ensembl.parse_line,
            'r': genes.refGene.parse_line,
            'k': genes.knownGene.parse_line,
            'a': genes.aceview.parse_line,
            'x': genes.ensg.parse_line
          }[model](line)
          if txt.name == txt_name or txt.alias == txt_name:
            txts.append(txt)
    return txts

  def SortedExons(self): #{
    return sorted((CoordPairCls(exon) for exon in self.exons),
      key=lambda exon: exon.min)
  #} end def
