#! /usr/bin/env python
"""
coord_pair.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import custom modules
from utils.error import MyError
from utils.messages import ErrMsg

class CoordPairCls: #{
  def __init__(self, *args, **kwargs): #{
    pos_strand = kwargs.get("pos_strand")
    if (0 == len(args)): #{
      self.start = None
      self.end   = None
      self.min = None
      self.max = None
      if (None == pos_strand):
        self.pos_strand = True
      else:
        self.pos_strand = pos_strand
      #} end if
    else:
      if (1 == len(args)): #{
        in_coords = args[0]
      else:
        in_coords = args
      #} end if
      if (isinstance(in_coords, str)): #{
        in_coords = in_coords.split("-")
      elif (isinstance(in_coords, int)):
        in_coords = (in_coords, in_coords)
      #} end if
      if (isinstance(in_coords, tuple) or
          isinstance(in_coords, list)):
        try: #{
          (self.start, self.end) = map(int, in_coords)
          sorted_coords = sorted((self.start, self.end))
          (self.min, self.max) = sorted_coords
        except ValueError, e:
          raise CoordPairError("cannot initialize coord pair from %s (%s)" %
            (in_coords, e))
        #} end try
      elif (isinstance(in_coords, CoordPairCls)):
        self.start = in_coords.start
        self.end   = in_coords.end
        self.min = in_coords.min
        self.max = in_coords.max
      else:
        raise CoordPairError("cannot initialize coord pair from %s (%s)" %
          (in_coords, type(in_coords)))
      #} end if
      self.pos_strand = (self.start < self.end)
      if (None != pos_strand and pos_strand != self.pos_strand): #{
        raise CoordPairError("Strand flag provided is inconsistent with "
          "coord order: pos_strand = %s, coords = %s-%s" % (pos_strand,
          self.start, self.end))
      #} end if
    #} end if
  #} end def

  def __len__(self): #{
    return self.Span()
  #} end def

  def __str__(self): #{
    return self.ToString()
  #} end def

  def __getitem__(self, key): #{
    if (0 == key): #{
      return self.start
    elif (1 == key):
      return self.end
    else:
      raise IndexError
    #} end if
  #} end def

  def Span(self): #{
    return (self.max-self.min)+1
  #} end def

  def Union(self, other): #{
    if (not isinstance(other, CoordPairCls)): #{
      other = CoordPairCls(other)
    #} end if
    if (None == self.min or other.min < self.min): #{
      self.min = other.min
    #} end if
    if (None == self.max or other.max > self.max): #{
      self.max = other.max
    #} end if
    if (self.pos_strand): #{
      (self.start, self.end) = (self.min, self.max)
    else:
      (self.start, self.end) = (self.max, self.min)
    #} end if
    return None
  #} end def

  def Contains(self, other): #{
    if (not isinstance(other, CoordPairCls)): #{
      other = CoordPairCls(other)
    #} end if
    if (self.min < other.min and other.max < self.max): #{
      return True
    #} end if
    return False
  #} end def

  def Overlaps(self, other): #{
    if (not isinstance(other, CoordPairCls)): #{
      other = CoordPairCls(other)
    #} end if
    if (other.max < self.min or self.max < other.min): #{
      return False
    #} end if
    return True
  #} end def

  def OverlapAmount(self, other): #{
    if (not isinstance(other, CoordPairCls)): #{
      other = CoordPairCls(other)
    #} end if
    #ErrMsg("BEFORE: %s" % self)
    intersect = self.Copy()
    try: #{
      intersect.Intersect(other)
      #ErrMsg("AFTER: %s" % self)
    except CoordPairError:
      #ErrMsg("AFTER: %s" % self)
      return 0
    #} end try
    return intersect.Span()
  #} end def

  def Intersect(self, other): #{
    if (not isinstance(other, CoordPairCls)): #{
      other = CoordPairCls(other)
    #} end if
    #ErrMsg("Intersecting %s with %s..." % (self, other))
    #ErrMsg("  Min check: %s" % (None == self.min or other.min > self.min))
    if (None == self.min or other.min > self.min): #{
      #ErrMsg("    UPDATE MIN!")
      self.min = other.min
    #} end if
    #ErrMsg("  Max check: %s" % (None == self.max or other.max < self.max))
    if (None == self.max or other.max < self.max): #{
      #ErrMsg("    UPDATE MAX!")
      self.max = other.max
    #} end if
    #ErrMsg("  After intersection: %s" % self)
    if (self.max < self.min): #{
      #ErrMsg("    NO INTERSECTION!")
      raise CoordPairError("coord-pairs do not intersect")
    #} end if
    if (self.pos_strand): #{
      (self.start, self.end) = (self.min, self.max)
    else:
      (self.start, self.end) = (self.max, self.min)
    #} end if
  #} end def

  def ToString(self): #{
    return "%i-%i" % (self.start, self.end)
  #} end def

  def Copy(self): #{
    return CoordPairCls(self)
  #} end def

  def ResortCoords(self): #{
    unsorted = (self.min, self.max)
    (self.min, self.max) = sorted(unsorted)
  #} end def

  def SetMin(self, new_min): #{
    self.min = new_min
    if (self.pos_strand): #{
      self.start = new_min
    else:
      self.end = new_min
    #} end if
  #} end def

  def SetMax(self, new_max): #{
    self.max = new_max
    if (self.pos_strand): #{
      self.end = new_max
    else:
      self.start = new_max
    #} end if
  #} end def

  def MoveMin(self, delta): #{
    self.min += delta
    if (self.pos_strand): #{
      self.start = self.min
    else:
      self.end = self.min
    #} end if
  #} end def

  def MoveMax(self, delta): #{
    self.max += delta
    if (self.pos_strand): #{
      self.end = self.max
    else:
      self.start = self.max
    #} end if
  #} end def
#} end class

def BetweenCoords(coords1, coords2): #{
  start = max(coords1.min, coords2.min)-1
  end = min(coords1.max, coords2.max)+1
  result = CoordPairCls((start, end))
  return result
#} end def

def CutOrExtendBlocks(input_blocks, target_length, from_left): #{
  output_blocks = list()
  remaining = target_length
  #ErrMsg("From left: %s. Initially, remaining = %i" % (from_left, remaining))
  if (not from_left): #{
    #ErrMsg("  (reversing block order)")
    input_blocks = reversed(input_blocks)
  #} end if
  for inblock in input_blocks: #{
    if (0 == remaining): #{
      break
    #} end if
    if (not isinstance(inblock, CoordPairCls)): #{
      inblock = CoordPairCls(inblock)
    #} end if
    #ErrMsg("Block: %s (%i)" % (inblock, inblock.Span()))
    if (remaining > inblock.Span()): #{
      output_blocks.append(inblock)
      remaining -= inblock.Span()
    else:
      if (from_left): #{
        part_block = CoordPairCls(inblock.min, (inblock.min+remaining)-1)
      else:
        part_block = CoordPairCls((inblock.max-remaining)+1, inblock.max)
      #} end if
      if (part_block.Span() != remaining): #{
        #ErrMsg("  part_block: %s (%i)" % (part_block, part_block.Span()))
        raise CoordPairError("improper calculation of partial block")
      #} end if
      output_blocks.append(part_block)
      remaining -= part_block.Span()
    #} end if
    #ErrMsg("  remaining = %i" % remaining)
  #} end for
  if (0 < remaining): #{
    last_block = output_blocks[-1]
    if (from_left): #{
      #last_block.max += remaining
      last_block.MoveMax(remaining)
    else:
      #last_block.min -= remaining
      last_block.MoveMin(-remaining)
    #} end if
  #} end if
  if (not from_left): #{
    return list(reversed(output_blocks))
  #} end if
  return output_blocks
#} end def

def MergeAdjacentBlocks(input_blocks): #{
  output_blocks = list()
  prev_block = None
  for inblock in input_blocks: #{
    if (not isinstance(inblock, CoordPairCls)): #{
      inblock = CoordPairCls(inblock)
    #} end if
    if (None != prev_block and prev_block.min > inblock.min): #{
      raise CoordPairCls("adjacent blocks cannot be merged when the blocks "
        "are out of order: %s" % input_blocks)
    #} end if
    if (0 == len(output_blocks) or
        output_blocks[-1].max < (inblock.min-1)): #{
      output_blocks.append(inblock)
    elif (output_blocks[-1].max == (inblock.min-1)):
      output_blocks[-1].SetMax(inblock.max)
    else:
      raise CoordPairCls("adjacent blocks cannot be merged when the blocks "
        "overlap: %s, %s" % (output_blocks[-1], inblock))
    #} end if
    prev_block = inblock
  #} end for
  return output_blocks
#} end def

#### EXCEPTION CLASSES ####
class CoordPairError(MyError): #{
  pass
#} end class
