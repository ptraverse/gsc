#! /usr/bin/env python
"""
gene_overlap.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import re

# import custom modules
from utils.error import MyError
from utils.general import RemoveIndexingFlag
from utils.messages import DebugMsg, ExtremeDebugMsg
from utils.multi_dict import SetupMultiDict
from parsers.tokenizer import TokenizerCls

# constants
MIN_FEATURE_OVERLAP = 5

class FeatureOverlapCls: #{
  def __init__(self, tokenizer, multi_target=False): #{
    #DebugMsg(tokenizer, "Parsing feature overlap line...")
    #DebugMsg(tokenizer, "Parsing chromosome")
    self.chrom = tokenizer.next()
    #DebugMsg(tokenizer, "Parsing start, end, and query_id")
    #(self.start, self.end, self.query_id) = map(int, tokenizer.next())
    for attr in ["start", "end"]: #{
      setattr(self, attr, int(tokenizer.next()))
    #} end for
    self.query_id = tokenizer.next()
    if (multi_target): #{
      self.ParseMultipleTargets(tokenizer)
    else:
      self.id_dict = dict()
      # overlaps have 6 columns, nearby regions have 8
      if (6 == tokenizer.Size()): #{
        self.ParseSingleTarget(tokenizer)
      elif (8 == tokenizer.Size()):
        self.ParseNearbyTarget(tokenizer)
      else:
        raise GeneFeatureOverlapError("wrong number of columns in overlap "
          "line: %s" % tokenizer.InputString())
      #} end if
    #} end if
  #} end def

  def ParseSingleTarget(self, tokenizer): #{
    #DebugMsg(tokenizer, "Parsing overlap results")
    ParseFeatureDescriptor(self.id_dict, tokenizer.next())
    self.amount = int(tokenizer.next())
    self.overlap_type = "overlap"
  #} end def

  def ParseNearbyTarget(self, tokenizer): #{
    #DebugMsg(tokenizer, "Parsing nearby results")
    #DebugMsg(tokenizer, "Parsing nearby pos")
    near_pos = int(tokenizer.next())
    #DebugMsg(tokenizer, "Parsing nearby descriptor")
    descriptor = tokenizer.next()
    #DebugMsg(tokenizer, "Parsing nearby distance")
    self.amount = int(tokenizer.next())
    #DebugMsg(self, "  Results: pos=%i,desc=%s,amount=%i" %
    #  (near_pos, descriptor, self.amount))
    feature_size = int(tokenizer.next())
    # check if no overlap or nearby gene was found
    if (0 == near_pos and "0" == descriptor and -1 == feature_size): #{
      self.overlap_type = "none"
      return
    #} end if
    ParseFeatureDescriptor(self.id_dict, descriptor)
    self.overlap_type = "nearby"
    #DebugMsg(self, "Finished parsing nearby results")
  #} end def

  def ParseMultipleTargets(self, tokenizer): #{
    self.target_list = list()
    targets = TokenizerCls(tokenizer.next(), delimiter="^")
    for target in targets: #{
      id_dict = dict()
      ParseFeatureDescriptor(id_dict, target)
      self.target_list.append(id_dict)
    #} end for
  #} end def

  def NumID(self): #{
    return int(self.query_id)
  #} end def
#} end class

def ParseFeatureDescriptor(id_dict, descriptor): #{
  descriptor_list = RemoveIndexingFlag(descriptor).split("~")
  (full_name, id_dict['feature'], id_dict['num']) = descriptor_list
  if (">" in full_name): #{
    (id_dict['transcript'], id_dict['gene']) = full_name.rsplit(">",1)
  else:
    id_dict['transcript'] = full_name
    id_dict['gene'] = "no_alias"
  #} end if
#} end def

class TranscriptOverlapCls: #{
  def __init__(self, id): #{
    self.id = id
    self.gene_name = None
    # features[feature_type] = feature_number_set
    self.features = dict()
    self.low_overlaps = dict()
    self.amount = 0
  #} end def

  def Update(self, overlap): #{
    if (self.id != overlap.id_dict['transcript']): #{
      raise GeneFeatureOverlapError("Attempting to update transcript overlap "
        "with feature with mismatching id: %s != %s" % (self.id,
        overlap.id_dict['transcript']))
    #} end if
    if (None == self.gene_name): #{
      self.gene_name = overlap.id_dict['gene']
    elif (self.gene_name != overlap.id_dict['gene']):
      raise GeneFeatureOverlapError("Mismatching gene names: %s != %s" %
        (self.gene_name, overlap.id_dict['gene']))
    #} end if
    if ("nearby" == overlap.overlap_type or
        overlap.amount > MIN_FEATURE_OVERLAP): #{
      if (overlap.id_dict['feature'] not in self.features): #{
        self.features[overlap.id_dict['feature']] = set()
      #} end if
      self.features[overlap.id_dict['feature']].update(overlap.id_dict['num'])
    else:
      if (overlap.id_dict['feature'] not in self.low_overlaps): #{
        self.low_overlaps[overlap.id_dict['feature']] = set()
      #} end if
      low_overlap_feature = self.low_overlaps[overlap.id_dict['feature']]
      low_overlap_feature.update(overlap.id_dict['num'])
    #} end if
    if ("overlap" == overlap.overlap_type): #{
      self.amount += overlap.amount
    elif ("nearby" == overlap.overlap_type):
      if (0 == self.amount): #{
        self.amount = overlap.amount
      else:
        self.amount = min(self.amount, overlap.amount)
      #} end if
    else:
      raise GeneFeatureOverlapError("Unrecognized overlap type: %s" %
        overlap.overlap_type)
    #} end if
  #} end def

  def CheckFeatures(self): #{
    if (0 == len(self.features)): #{
      self.features = self.low_overlaps
    else:
      for feature in self.features: #{
        if (0 == len(self.features[feature]) and
            feature in self.low_overlaps): #{
          self.features[feature] = self.low_overlaps[feature]
        #} end if
      #} end for
  #} end def
#} end class

def AddOverlapToAlign(align, overlap): #{
  if (not hasattr(align, "overlap")): #{
    align.overlap = dict()
    align.best_transcripts_chosen = False
  #} end if
  if (not hasattr(align, "nearby")): #{
    align.nearby = dict()
  #} end if
  if ("none" == overlap.overlap_type): #{
    return
  #} end if
  overlap_dict = getattr(align, overlap.overlap_type)
  if (overlap.id_dict['transcript'] not in overlap_dict): #{
    new_transcript = TranscriptOverlapCls(overlap.id_dict['transcript'])
    overlap_dict[overlap.id_dict['transcript']] = new_transcript
  #} end if
  overlap_dict[overlap.id_dict['transcript']].Update(overlap)
  if (not hasattr(align, "best_"+overlap.overlap_type)): #{
    best = None
  else:
    best = getattr(align, "best_"+overlap.overlap_type)
  #} end if
  new_amount = overlap_dict[overlap.id_dict['transcript']].amount
  if (None == best or best < new_amount): #{
    setattr(align, "best_"+overlap.overlap_type, new_amount)
  #} end if
#} end def

def ChooseBestTranscripts(align, buffer): #{
  # if the best transcripts have already been chosen, do not do anything
  if (align.best_transcripts_chosen): #{
    return
  #} end if
  # if there are any true overlaps, discard nearby genes
  if (hasattr(align, "overlap") and None != align.overlap and
      0 < len(align.overlap)): #{
    align.nearby = dict()
    good_transcripts = list()
    for transcript in align.overlap.itervalues(): #{
      if ((align.best_overlap - buffer) < transcript.amount): #{
        good_transcripts.append(transcript)
      #} end if
    #} end for
    # overlap[gene][transcript][feature] = feature number set
    align.overlap = dict()
    for transcript in good_transcripts: #{
      transcript.CheckFeatures()
      keys = [transcript.gene_name, transcript.id]
      SetupMultiDict(align.overlap, keys, transcript.features, list_end=False)
    #} end for
  else:
    align.overlap = dict()
    if (hasattr(align, "nearby") and None != align.nearby and
        0 < len(align.nearby)): #{
      good_transcripts = list()
      for transcript in align.nearby.itervalues(): #{
        if ((align.best_nearby + buffer) > transcript.amount): #{
          good_transcripts.append(transcript)
        #} end if
      #} end for
      # nearby[gene][transcript] = distance
      align.nearby = dict()
      for transcript in good_transcripts: #{
        keys = [transcript.gene_name, transcript.id]
        SetupMultiDict(align.nearby, keys, transcript.amount, list_end=False)
      #} end for
    else:
      align.nearby = dict()
    #} end if
  #} end if
  align.best_transcripts_chosen = True
#} end def

def GenesOverlapped(align, add_gene_annotation): #{
  if (add_gene_annotation): #{
    if (not hasattr(align, "overlap") or None == align.overlap): #{
      raise AlignmentError("Cannot get gene features overlapped for "
        "alignment without gene feature overlap information")
    #} end if
    return GeneOverlapString(align.overlap)
  else:
    return "N/A"
  #} end if
#} end def

def GeneOverlapString(gene_dict): #{
  if (0 < len(gene_dict)): #{
    gene_strs = list()
    for gene_name in sorted(gene_dict.keys(), key=str.lower): #{
      #if ("no_alias" == gene_name): #{
      #  continue
      #} end if
      transcript_dict = gene_dict[gene_name]
      gene_parts = set()
      non_coding = True
      for transcript_id in transcript_dict.iterkeys(): #{
        gene_parts.update(transcript_dict[transcript_id].keys())
        if (not transcript_id.endswith("non_coding")): #{
          non_coding = False
        #} end if
      #} end for
      if (non_coding): #{
        gene_name += ".non_coding"
      #} end if
      #transcript_str = TranscriptOverlapString(gene_dict[gene_name])
      parts_str = ",".join(sorted(gene_parts, key=str.lower))
      #gene_strs.append("%s(%s)" % (gene_name, transcript_str))
      gene_strs.append("%s(%s)" % (gene_name, parts_str))
    #} end for
    #if ("no_alias" in gene_dict): #{
    #  transcript_str = TranscriptOverlapString(gene_dict["no_alias"])
    #  gene_strs.append("%s" % (transcript_str))
    #} end for
    return ",".join(gene_strs)
  else:
    return "none"
  #} end if
#} end def

def TranscriptOverlapString(transcript_dict): #{
  transcript_strs = list()
  for transcript_id in sorted(transcript_dict.keys(), key=str.lower): #{
    #feature_str = FeatureOverlapString(transcript_dict[transcript_id])
    feature_str = ",".join(sorted(transcript_dict[transcript_id],
      key=str.lower))
    transcript_strs.append("%s(%s)" % (transcript_id, feature_str))
  #} end for
  return ",".join(transcript_strs)
#} end def

def FeatureOverlapString(feature_dict): #{
  feature_strs = list()
  for feature_type in sorted(feature_dict.keys(), key=str.lower): #{
    nums_str = ",".join(sorted(feature_dict[feature_type], key=int))
    feature_strs.append("%s(%s)" % (feature_type, nums_str))
  #} end for
  return ",".join(feature_strs)
#} end def

def NearbyGenes(align, add_gene_annotation): #{
  if (add_gene_annotation): #{
    if (not hasattr(align, "nearby") or None == align.nearby): #{
      raise AlignmentError("Cannot get nearby genes for alignment without "
        "gene feature overlap information")
    #} end if
    return GeneNearbyString(align.nearby)
  else:
    return "N/A"
  #} end if
#} end def

def GeneNearbyString(gene_dict): #{
  if (0 < len(gene_dict)): #{
    gene_strs = list()
    for gene_name in sorted(gene_dict.keys(), key=str.lower): #{
      #if ("no_alias" == gene_name): #{
      #  continue
      #} end if
      transcript_strs = list()
      for transcript_id in sorted(gene_dict[gene_name], key=str.lower): #{
        distance =  gene_dict[gene_name][transcript_id]
        transcript_strs.append("%s(%i)" % (transcript_id, distance))
      #} end for
      gene_strs.append("%s(%s)" % (gene_name, ",".join(transcript_strs)))
    #} end for
    #if ("no_alias" in gene_dict): #{
    #  transcript_strs = list()
    #  for transcript_id in sorted(gene_dict[gene_name], key=str.lower): #{
    #    distance =  gene_dict[gene_name][transcript_id]
    #    transcript_strs.append("%s(%i)" % (transcript_id, distance))
    #  #} end for
    #  gene_strs.append("%s(%s)" % (gene_name, ",".join(transcript_strs)))
    #} end for
    return ",".join(gene_strs)
  else:
    return "none"
  #} end if
#} end def

#### EXCEPTION CLASSES ####
class GeneFeatureOverlapError(MyError): #{
  pass
#} end class
