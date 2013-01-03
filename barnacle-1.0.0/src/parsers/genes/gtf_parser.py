#! /usr/bin/env python
"""
gtf_parser.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
from utils.messages import DebugMsg, ExtremeDebugMsg
from utils.files_paths import FileBoxCls
from parsers.tokenizer import TokenizerCls
from transcript import Transcript

# CONSTANTS

class GTFAnnotationParserCls: #{
  def __init__(self, input_path, log_info=None): #{
    self.file = FileBoxCls(input_path, "r",
      "cannot read gene annotations input file")
    self.curr_feature = None
    self.log_info = log_info
    self.finished = False
  #} end def

  def __del__(self): #{
    self.close()
  #} end def

  def __iter__(self): #{
    return self
  #} end def

  def next(self): #{
    if (self.finished): #{
      raise StopIteration
    #} end if
    transcript = None
    try: #{
      if (None == self.curr_feature): #{
        self.ParseFeature()
      #} end if
      transcript = GTFTranscriptCls(name=self.curr_feature.name)
      while (self.curr_feature.name == transcript.name): #{
        transcript.Update(self.curr_feature)
        self.ParseFeature()
      #} end while
    except StopIteration:
      self.finished = True
    #} end try
    if (None == transcript): #{
      raise StopIteration
    #} end if
    transcript.CreateExonList()
    return transcript
  #} end def

  def ParseFeature(self): #{
    #ExtremeDebugMsg(self, "  Parsing feature from file...")
    try: #[
      line = self.file.next()
    except StopIteration, e:
      self.curr_feature = None
      raise e
    #} end try
    tokenizer = TokenizerCls(line, delimiter="\t", log_info=self.log_info)
    self.curr_feature = GTFFeatureCls(tokenizer)
  #} end def

  def Reset(self): #{
    self.close()
    self.curr_feature = None
  #} end def

  def close(self): #{
    if (hasattr(self, "file") and None != self.file): #{
      self.file.close()
    #} end if
  #} end def
#} end class

class GTFFeatureCls: #{
  def __init__(self, tokenizer): #{
    self.chrom  = tokenizer.next()
    self.t_type = tokenizer.next()
    self.type   = tokenizer.next()
    self.start  = int(tokenizer.next())
    self.end    = int(tokenizer.next())
    # skip the score
    tokenizer.Skip()
    self.strand = tokenizer.next()
    # skip the frame
    tokenizer.Skip()
    # self.name = transcript ID
    self.name     = None
    self.exon_num = None
    # self.alias = gene name
    self.alias    = None
    # retokenize the remaining attributes
    attributes = TokenizerCls(tokenizer.next(), delimiter=";",
      log_info=tokenizer.log_info)
    for attribute in attributes: #{
      #DebugMsg(tokenizer, "Attribute: \"%s\"" % attribute)
      try: #{
        (type, value) = attribute.strip(" ").split(" ", 1)
      except ValueError, e:
        raise GTFAnnotationParserError("Cannot parse attribute %s: %s" %
          (attribute, e))
      #} end try
      value = value.strip("\"")
      if ("transcript_id" == type): #{
        self.name = value
      elif ("exon_number" == type):
        self.exon_num = int(value)-1
      elif ("gene_name" == type):
        self.alias = value
      elif (type in ["gene_id","transcript_name","protein_id"]):
        # ignore these attributes
        continue
      # add code to handle "gene_biotype" field
      else:
        DebugMsg(tokenizer, "Unrecognized attribute type \"%s: %s\"" %
          (type, value))
      #} end if
    #} end for
  #} end def
#} end class

class GTFTranscriptCls(Transcript): #{
  def __init__(self, name=None): #{
    self.chrom    = None
    self.t_type   = None
    self.exon_d   = dict()
    self.exons    = list()
    self.txStart  = None
    self.txEnd    = None
    self.cdsStart = None
    self.cdsEnd   = None
    # self.name = transcript ID
    self.name     = name
    # self.alias = gene name
    self.alias    = None
    self.strand   = None
  #} end def

  def Update(self, feature): #{
    errors = list()
    if (None == self.chrom): #{
      self.chrom = feature.chrom
    elif (feature.chrom != self.chrom):
      errors.append("chroms do not match: \"%s\" - \"%s\"" %
        (self.chrom, feature.chrom))
    #} end if
    if (None == self.t_type): #{
      self.t_type = feature.t_type
    elif (feature.t_type != self.t_type):
      errors.append("t_types do not match: \"%s\" - \"%s\"" %
        (self.t_type, feature.t_type))
    #} end if
    if (None == self.name): #{
      self.name = feature.name
    elif (feature.name != self.name):
      errors.append("names do not match: \"%s\" - \"%s\"" %
        (self.name, feature.name))
    #} end if
    if (None == self.alias): #{
      self.alias = feature.alias
    elif (feature.alias != self.alias):
      errors.append("aliases do not match: \"%s\" - \"%s\"" %
        (self.alias, feature.alias))
    #} end if
    if (None == self.strand): #{
      self.strand = feature.strand
    elif (feature.strand != self.strand):
      errors.append("strands do not match: \"%s\" - \"%s\"" %
        (self.strand, feature.strand))
    #} end if
    if ("exon" == feature.type): #{
      f_coords = (feature.start, feature.end)
      if (feature.exon_num in self.exon_d): #{
        e_coords = self.exon_d[feature.exon_num]
        if (f_coords != e_coords): #{
          raise GTFAnnotationParserError("conflicting exon coords for %s exon "
            "%i: %s-%s, %s-%s" % (self.name, feature.exon_num, e_coords[0],
            e_coords[1], f_coords[0], f_coords[1]))
        #} end if
      else:
        self.exon_d[feature.exon_num] = f_coords
        if (None == self.txStart or feature.start < self.txStart): #{
          self.txStart = feature.start
        #} end if
        if (None == self.txEnd or feature.end > self.txEnd): #{
          self.txEnd = feature.end
        #} end if
      #} end if
    elif ("CDS" == feature.type):
      if (None == self.cdsStart or feature.start < self.cdsStart): #{
        self.cdsStart = feature.start
      #} end if
      if (None == self.cdsEnd or feature.end > self.cdsEnd): #{
        self.cdsEnd = feature.end
      #} end if
    # ignore these features
    elif (feature.type not in ["start_codon", "stop_codon"]):
      raise GTFAnnotationParserError("Unrecognized feature type for gene %s: "
        "%s" % (feature.name, feature.type))
    #} end if
    if (0 < len(errors)): #{
      raise GTFAnnotationParserError("Error updating gene %s\n  %s" %
        (self.alias, "\n  ".join(errors)))
    #} end if
  #} end def

  def CreateExonList(self): #{
    if (None == self.cdsStart): #{
      self.cdsStart = self.txStart
    #} end if
    if (None == self.cdsEnd): #{
      self.cdsEnd = self.txStart
    #} end if
    for i in range(max(self.exon_d.keys())+1): #{
      if (i in self.exon_d): #{
        self.exons.append(self.exon_d[i])
      else:
        self.exons.append(None)
      #} end if
    #} end for
  #} end def

  #def FullName(self): #{
  #  if (None == self.full_name): #{
  #    self.full_name = ".".join(self.t_id, self.gene_name)
  #  #} end if
  #  return self.full_name
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class GTFAnnotationParserError(MyError): #{
  pass
#} end class
