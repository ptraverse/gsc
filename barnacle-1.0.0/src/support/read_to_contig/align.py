#! /usr/bin/env python
"""
align.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import re

# import custom modules
from utils.error import MyError
from utils.messages import LogMsg, TestMsg
from support.SAM_record import SAM_record

class R2CAlignCls: #{
  def __init__(self, source, log_info=None): #{
    self.id            = None # qname
    self.ctg_id        = None
    self.length        = None # qlen
    self.left          = None # pos + 1 (convert from 0-based to 1-based)
    self.right         = None # left + length
    self.strand        = None
    self.num_best_hits = None # opt('X0')
    self.num_subopt    = 0    # opt('X1'), presume 0 because sometimes not present in record
    self.cigar         = None # convert cigar to string
    self.edit_dist     = None # opt('NM')
    self.is_perfect    = None
    self.type          = None # opt('XT') -- Unique/Repeat/N/Mate-sw
    self.alt_hits      = list() # opt('XA') -- chr,pos,CIGAR,NM;
    self.curr_alt      = -1
    self.alt_missing   = False # if param was too low to output alternate hits
    self.is_alt_hit    = False
    self.opt           = dict()
    self.log_info      = log_info
    if (isinstance(source, SAM_record)): #{
      self.InitializeFromSAMRecord(source)
    elif (isinstance(source, R2CAlignCls)):
      self.InitializeFromParentAlign(source)
    else:
      raise R2CAlignError("Unrecognized read alignment source type: %s" %
        source)
    #} end if
  #} end def

  # delay determining "unique"
  #def InitializeFromSAMRecord(self, record, ctgs_in_event): #{
  def InitializeFromSAMRecord(self, record): #{
    record.sanitize_qname()
    self.ctg_id        = CleanAlignContigID(record.rname)
    self.id            = record.qname
    if (None == record.read): #{
      LogMsg(self.log_file,
        "Warning read %s is neither first nor second read" % self.id)
    else:
      self.id += "_%i" % record.read
    #} end if
    # try to get length from CIGAR
    cigar_match = re.search(r"([0-9]+S)?(?P<M>[0-9]+)M([0-9]+S)?",
        record.cigar)
    if (None == cigar_match): #{
      self.length = record.read_length
    else:
      self.length = int(cigar_match.group("M"))
    #} end if
    if (None == self.length or 0 == self.length): #{
      raise R2CAlignError("Could not get read length for read %s" % self.id)
    #} end if
    self.left          = record.pos # already 1-based
    self.right         = self.left + self.length - 1 # assume no gaps
    if ("Reverse" == record.strand): #{
      self.strand = "-"
    elif ("Forward" == record.strand):
      self.strand = "+"
    else:
      raise R2CAlignError("Unrecognized strand \"%s\" for read %s" %
        (record.strand, self.id))
    #} end if
    self.GetOptFields(record.opt_fields)
    self.num_best_hits = self.opt.get('X0',1)
    self.num_subopt    = self.opt.get('X1',0)
    self.cigar         = record.cigar
    self.edit_dist     = record.edit_distance
    self.is_perfect    = record.perfect_match()
    self.type          = self.opt.get('XT')
    if ('XA' in self.opt): #{
      self.alt_hits.extend(self.opt.get('XA').rstrip(";").split(";"))
    elif (1 < self.num_best_hits):
      self.alt_missing = True
    #} end if
  #} end def

  def GetOptFields(self, opt_field_list): #{
    for opt_field_str in opt_field_list: #{
      try: #{
        (tag, type, value) = opt_field_str.split(":", 2)
      except ValueError, e:
        raise R2CAlignError("Cannot parse optional field: \"%s\"\n%s" %
          (opt_field_str, e))
      #} end try
      if ("XA" == tag and type not in ["A","Z"]): #{
        raise R2CAlignError("Invalid field type for alternative hits field "
          "(XA): \"%s\", should be \"A\" or \"Z\"." % type)
      #} end if
      if (type in ["A","Z"]): #{
        self.opt[tag] = value
      elif ("i" == type):
        self.opt[tag] = int(value)
      elif ("f" == type):
        self.opt[tag] = float(value)
      else:
        raise R2CAlignError("Unrecognized optional field type "
          "\"%s\" for read %s" % (type, self.id))
      #} end if
    #} end for
  #} end def

  def InitializeFromParentAlign(self, parent): #{
    alt_hit = parent.alt_hits[parent.curr_alt]
    (ctg_id, full_pos, self.cigar, edit_dist) = alt_hit.split(",")
    self.ctg_id        = CleanAlignContigID(ctg_id)
    self.id            = parent.id
    self.length        = parent.length
    self.left          = int(full_pos[1:]) # already 1-based
    self.right         = self.left + self.length - 1
    self.strand        = full_pos[0]
    self.num_best_hits = parent.num_best_hits
    self.num_subopt    = parent.num_subopt
    self.edit_dist     = int(edit_dist)
    # a perfect match has every base as M in the CIGAR string and
    #   an edit distance of 0 (since M can be a mismatch)
    self.is_perfect    = ("%iM" % self.length == self.cigar and
      0 == self.edit_dist)
    self.type          = parent.type
    self.is_alt_hit    = True
    self.alt_hits.extend(parent.alt_hits[:parent.curr_alt])
    self.alt_hits.extend(parent.alt_hits[parent.curr_alt+1:])
    self.alt_hits.append(parent.AltString())
  #} end def

  def IsUnique(self, ctgs_in_event): #{
    TestMsg("Testing whether read alignment is unique. Num best hits: %i" %
        self.num_best_hits + " Contigs in event: %i" % len(ctgs_in_event))
    if (1 == self.num_best_hits): #{
      TestMsg("UNIQUE: 1 best hit.")
      return True
    elif (self.num_best_hits < (2 * len(ctgs_in_event))+1):
      alt_hit_ctgs = set()
      for alternate_hit in self.alt_hits: #{
        raw_ctg_id = alternate_hit.split(",")[0]
        ctg_id = CleanAlignContigID(raw_ctg_id)
        alt_hit_ctgs.add(ctg_id)
      #} end for
      # add the primary hit contig id
      raw_ctg_id = self.ctg_id
      ctg_id = CleanAlignContigID(raw_ctg_id)
      alt_hit_ctgs.add(ctg_id)
      alt_hit_ctgs_in_event = alt_hit_ctgs.intersection(ctgs_in_event)
      fraction_in_event = (float(len(alt_hit_ctgs_in_event)) /
        float(len(alt_hit_ctgs)))
      TestMsg("\n".join([
        "Alt hit contigs:  %s" % ",".join(sorted(alt_hit_ctgs)),
        "Contigs in event: %s" % ",".join(sorted(ctgs_in_event)),
        "Alt in event:     %s" % ",".join(sorted(alt_hit_ctgs_in_event)),
        "Fraction: %.3f" % fraction_in_event]))
      if (0.5 <= fraction_in_event): #{
        TestMsg("UNIQUE: hits mostly in event")
        return True
      #} end if
    #} end if
    TestMsg("NOT unique: too many best hits.")
    return False
  #} end def

  def ToString(self): #{
    if (None == self.edit_dist): #{
      edit_dist = -1
    else:
      edit_dist = self.edit_dist
    #} end if
    data = [
      "ID:%s" % self.id,
      "CTG:%s" % self.ctg_id,
      "COORDS:%i-%i(%ibp)" % (self.left, self.right, self.length),
      "CIGAR:%s" % self.cigar,
      "EDIT_DIST:%i" % edit_dist,
      "PERFECT:%s" % self.is_perfect,
      "TYPE:%s" % self.type,
      "NUM_HITS:%i" % self.num_best_hits,
      "ALT_MISSING:%s" % self.alt_missing,
    ]
    if (0 == len(self.alt_hits)): #{
      if (self.alt_missing): #{
        data.append("ALT:missing")
      else:
        data.append("ALT:none")
      #} end if
    else:
      data.append("ALT:%s" % ";".join(self.alt_hits))
    #} end if
    return " ".join(data)
  #} end def

  def AltString(self): #{
    alt_data_list = [self.ctg_id, "%s%i" % (self.strand, self.left),
      self.cigar, "%i" % self.edit_dist]
    return ",".join(alt_data_list)
  #} end def
#} end class

def CleanAlignContigID(align_ctg_id):
  # remove library name, replace any symbols between k-value and
  #   numeric portion of id with "_"
  return re.sub(r"^[^k]*(k[0-9]+)[^0-9]*", r"\1_", align_ctg_id)
# end def

#### EXCEPTION CLASSES ####
class R2CAlignError(MyError): #{
  pass
#} end class
