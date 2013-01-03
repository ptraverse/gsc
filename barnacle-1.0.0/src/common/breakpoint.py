#! /usr/bin/env python
"""
breakpoint.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
import re

# import custom modules
from utils.error import ParsingError
from utils.general import ChrToInt

class BreakpointCls: #{
  def __init__(self, breakpoint_str): #{
    breakpoint_pattern = (r"(chr)?(?P<chr>[0-9]+|[XYM]):(?P<coord>[0-9]+)"
      "\((?P<dir>up|down)\)")
    breakpoint_match = re.search(breakpoint_pattern, breakpoint_str)
    if (None == breakpoint_match): #{
      raise ParsingError("cannot parse breakpoint coordinate: %s" %
        breakpoint_str)
    #} end if
    self.chr   = breakpoint_match.group("chr")
    self.coord = int(breakpoint_match.group("coord"))
    self.dir   = breakpoint_match.group("dir")
  #} end def

  def __str__(self): #{
    return self.ToString()
  #} end def

  def GetPositionKey(self, buffer): #{
    # pos_key = floor(pos / big_buffer)
    return (self.coord / buffer)
  #} end def

  def ToString(self): #{
    return "chr%s:%i(%s)" % (self.chr, self.coord, self.dir)
  #} end def
#} end def

def SortBreakpointTuple(breakpoints): #{
  return SortBreakpoints(breakpoints[0], breakpoints[1])
#} end def

def SortBreakpoints(breakpointA, breakpointB): #{
  if ((ChrToInt(breakpointA.chr) < ChrToInt(breakpointB.chr)) or
      (breakpointA.chr == breakpointB.chr and
       breakpointA.coord  < breakpointB.coord) or
      (breakpointA.chr == breakpointB.chr and
       breakpointA.coord == breakpointB.coord and
       "up" == breakpointA.dir)): #{
    return (breakpointA, breakpointB)
  else:
    return (breakpointB, breakpointA)
  #} end if
#} end def

def BreakpointsAreEqual(breakpointA, breakpointB, buffer=0,
    check_dirs=False): #{
  if (breakpointA.chr != breakpointB.chr): #{
    return False
  #} end if
  if (buffer < abs(breakpointA.coord - breakpointB.coord)): #{
    return False
  #} end if
  if (check_dirs and breakpointA.dir != breakpointB.dir): #{
    return False
  #} end if
  return True
#} end def

def BreakpointPairsAreEqual(pairA, pairB, buffer=0,
    check_dirs=False): #{
  pairA_sorted = SortBreakpoints(pairA[0], pairA[1])
  pairB_sorted = SortBreakpoints(pairB[0], pairB[1])
  if (not BreakpointsAreEqual(pairA_sorted[0], pairB_sorted[0],
      buffer, check_dirs)): #{
    return False
  #} end if
  if (not BreakpointsAreEqual(pairA_sorted[1], pairB_sorted[1],
      buffer, check_dirs)): #{
    return False
  #} end if
  return True
#} end def

def BreakpointPairString(pair, sort=False): #{
  if (sort): #{
    local_pair = SortBreakpoints(pair[0],pair[1])
  else:
    local_pair = pair
  #}
  return "%s-%s" % (local_pair[0].ToString(), local_pair[1].ToString())
#} end def
