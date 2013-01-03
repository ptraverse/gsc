#! /usr/bin/env python
"""
breakpoint_pair.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from utils.error import MyError

# import custom modules
from common.breakpoint import SortBreakpoints, BreakpointPairsAreEqual

# CONSTANTS

class BreakpointPairCls: #{
  def __init__(self, buffer): #{
    #print "Creating bp pair object"
    self.buffer = buffer
    self.pair = None
    self.barnacle_preds = list()
    self.other_preds = list()
    self.other_quality = None
  #} end def

  def AddBarnaclePrediction(self, b_member): #{
    new_pair = SortBreakpoints(b_member.breakpointA, b_member.breakpointB)
    #print "PAIR: ", self.pair
    if (None == self.pair): #{
      #print "Setting pair to %s-%s" % (new_pair[0].ToString(),
      #  new_pair[1].ToString())
      self.pair = new_pair
    elif (not BreakpointPairsAreEqual(self.pair, new_pair, self.buffer)): #{
      raise BreakpointPairError("Barnacle breakpoint does not match: %s, %s" %
        (BreakpointPairString(self.pair), BreakpointPairString(new_pair)))
    #} end if
    self.barnacle_preds.append(b_member.IDString())
  #} end def

  def AddOtherPrediction(self, pred): #{
    new_pair = SortBreakpoints(pred.breakpointA, pred.breakpointB)
    if (None == self.pair): #{
      self.pair = new_pair
    elif (not BreakpointPairsAreEqual(self.pair, new_pair, self.buffer)): #{
      raise BreakpointPairError("Barnacle breakpoint does not match: %s, %s" %
        (BreakpointPairString(self.pair), BreakpointPairString(new_pair)))
    #} end if
    self.other_preds.append(pred.id)
    if (hasattr(pred, "quality")): #{
      if (None == self.other_quality or
          self.other_quality < pred.quality): #{
        self.other_quality = pred.quality
      #} end if
    #} end if
  #} end def

  def Common(self): #{
    return (0 < len(self.barnacle_preds) and 0 < len(self.other_preds))
  #} end def

  def BarnacleOnly(self): #{
    return (0 < len(self.barnacle_preds) and 0 == len(self.other_preds))
  #} end def

  def OtherOnly(self): #{
    return (0 == len(self.barnacle_preds) and 0 < len(self.other_preds))
  #} end def

  def ToString(self): #{
    if (None == self.pair): #{
      pair_data = "Pair data not set!"
    else:
      pair_data = "%s-%s" % (self.pair[0].ToString(), self.pair[1].ToString())
    #} end if
    if (0 == len(self.barnacle_preds)): #{
      barnacle_str = "None"
    else:
      barnacle_str = ",".join(self.barnacle_preds)
    #} end if
    if (0 == len(self.other_preds)): #{
      other_str = "None"
    else:
      other_str = ",".join(self.other_preds)
    #} end if
    if (None == self.other_quality): #{
      qual_str = "None"
    else:
      qual_str = "%i" % (self.other_quality)
    #} end if
    return "\t".join([pair_data, barnacle_str, other_str, qual_str])
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class BreakpointPairError(MyError): #{
  pass
#} end class
