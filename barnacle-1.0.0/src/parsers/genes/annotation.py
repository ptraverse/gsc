#! /usr/bin/env python
"""
annotation.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules
from utils.error import MyError
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import FileBoxCls
from gtf_parser import GTFAnnotationParserCls
import aceview, ensembl, ensg, knownGene, refGene

# CONSTANTS
EXTS = {
  "gtf":"gtf",
}
ANNOT_TYPES = list([
  'aceview',
  'ensembl',
  'knowngene',
  'refseq',
])
PARSERS = {
  'gtf': GTFAnnotationParserCls,
}
PARSE_FUNCTIONS = {
  'aceview':   aceview.parse_line,
  'ensembl':   ensembl.parse_line,
  'ensg':      ensg.parse_line,
  'knowngene': knownGene.parse_line,
  'refseq':    refGene.parse_line,
}

class GeneAnnotationParserCls: #{
  def __init__(self, path, type=None, log_info=None): #{
    if (None == type): #{
      type = GetAnnotationsType(path)
    #} end if
    if (type in PARSERS): #{
      self.parser = PARSERS[type](path, log_info=log_info)
      self.file = None
      self.ParseLine = None
    elif (type in PARSE_FUNCTIONS): #{
      self.parser = None
      self.file = FileBoxCls(path, "r",
      "cannot open %s annotations file" % type)
      self.ParseLine = PARSE_FUNCTIONS[type]
    else:
      raise GeneAnnotationError("cannot determine correct annotation parser "
        "from annotations type: %s" % type)
    #} end if
    self.log_info = log_info
    self.finished = False
  #} end def

  def __del__(self): #{
    self.close()
  #} end def

  def __iter__(self): #{
    #if (None == self.parser): #{
    #  return self
    #else:
    #  return self.parser
    #} end if
    return self
  #} end def

  def next(self): #{
    if (self.finished): #{
      raise StopIteration
    #} end if
    #ExtremeDebugMsg(self, "Parsing annotation from file...")
    transcript = None
    try:
      if (None == self.parser): #{
        #ExtremeDebugMsg(self, "Using ParseLine function...")
        line = self.file.next()
        transcript = self.ParseLine(line)
      else:
        #ExtremeDebugMsg(self, "Using internal parser...")
        transcript = self.parser.next()
      #} end if
    except StopIteration:
      self.finished = True
    #} end try
    if (None == transcript): #{
      raise StopIteration
    #} end if
    transcript.gene_name = transcript.alias.replace(" ","_")
    transcript.transcript_id = transcript.name.replace(" ","_")
    #ExtremeDebugMsg(self, "Parsing transcript: %s (%s)" %
    #  (transcript.gene_name, transcript.transcript_id))
    return transcript
  #} end def

  def Close(self): #{
    for attr in ["file", "parser"]: #{
      if (hasattr(self, attr) and None != getattr(self, attr)): #{
        getattr(self, attr).close()
      #} end if
    #} end for
  #} end def

  def close(self): #{
    self.Close()
  #} end def
#} end class

def GetAnnotationsType(paths): #{
  # make sure paths is a list
  if (not isinstance(paths, list)): #{
    if (isinstance(paths, str)): #{
      if ("," in paths): #{
        paths = paths.split(",")
      elif (";" in paths):
        paths = paths.split(";")
      else:
        paths = [paths]
      #} end if
    elif (isinstance(paths, (set,tuple))):
      paths = list(paths)
    elif (isinstance(paths, dict)):
      paths = paths.values()
    elif (not hasattr(paths, '__iter__')):
      paths = [paths]
    else:
      raise GeneAnnotationError("cannot create list from paths object in "
         "GetAnnotationsType: %s (%s)" % (str(paths), type(paths).__name__))
    #} end if
  #} end if
  # check the paths the easy way (extensions first)
  for ext in EXTS: #{
    for path in paths: #{
      if (path.lower().endswith(ext.lower())): #{
        return EXTS[ext]
      #} end if
    #} end for
  #} end for
  # then full path names
  for annot_type in ANNOT_TYPES: #{
    for path in paths: #{
      if (annot_type.lower() in path.lower()): #{
        return annot_type
      #} end if
    #} end for
  #} end for
  # check for special file names
  for path in paths: #{
    if ("acembly" in path.lower()): #{
      return "aceview"
    elif ("ensgene" in path.lower()):
      return "ensembl"
    elif ("refgene" in path.lower()):
      return "refseq"
    elif ("ucsc_genes" in path.lower()):
      return "knowngene"
    #} end if
  #} end for
  # if you have not found an annotations type by this point, raise an error
  raise GeneAnnotationError("could not determine annotations type")
#} end def

#### EXCEPTION CLASSES ####
class GeneAnnotationError(MyError): #{
  pass
#} end class
