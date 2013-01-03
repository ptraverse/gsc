#! /usr/bin/env python
"""
group.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from operator import attrgetter

# import custom modules
from utils.error import MyError
from utils.general import GetGroupID
from utils.messages import ExtremeDebugMsg

class R2CGroupCls: #{
  def __init__(self, group_id, log_info=None): #{
    self.log_info = log_info
    self.group_id = int(group_id)
    self.members = dict()
    self.contigs = set()
  #} end def

  def AddMember(self, member_id, ctg_id, region_coords=None, mates=None): #{
    self.members[member_id] = R2CMemberCls(member_id, ctg_id,
      region_coords, log_info=self.log_info)
    self.contigs.add(ctg_id)
    self.AddMates(mates)
  #} end def

  def AddMates(self, mates): #{
    if (None == mates or "none" == mates.lower()): #{
      return
    #} end if
    for mate_str in mates.split(";"): #{
      (ctg_id, member_id, region_coords) = mate_str.split(",")
      if (member_id in self.members): #{
        raise R2CGroupError \
          ("duplicated member encountered: %s" % member_id)
      #} end if
      self.AddMember(member_id, ctg_id, region_coords)
    #} end for
  #} end def

  def InitializeWindows(self, min_overlap): #{
    for member in self.members.itervalues(): #{
      member.InitializeWindows(min_overlap)
    #} end for
  #} end def

  def AddSupport(self, from_member): #{
    if (from_member.member_id in self.members): #{
      self.members[from_member.member_id].AddSupport(from_member)
    else:
      self.members[from_member.member_id] = from_member
    #} end if
  #} end def

  def DebugString(self): #{
    data_list = list()
    sorted_members = \
      sorted(self.members.itervalues(), key=attrgetter('member_id'))
    for member_id in sorted(self.members.keys()): #{
      data_list.append(self.members[member_id].DebugString())
    #} end for
    return "\n".join(data_list)
  #} end def
#} end class

class R2CMemberCls: #{
  def __init__(self, member_id, ctg_id=None, region_coords=None,
      log_info=None): #{
    self.log_info = log_info
    self.ctg_id = ctg_id
    self.group_id = GetGroupID(member_id)
    self.member_id = member_id
    if (None == region_coords): #{
      self.left = self.right = None
    else:
      (self.left, self.right) = map(int, region_coords.split("-"))
      if (self.left > self.right): #{
        raise R2CGroupError("Invalid event region coordinates: %s" %
          region_coords)
      #} end if
    #} end if
    self.event_windows = None
  #} end def

  def Span(self):
    return (self.right - self.left) + 1
  # end def

  def InitializeWindows(self, min_overlap): #{
    self.event_windows = list()
    window_width = (2 * min_overlap) + 1
    if (self.Span() < window_width): #{
      new_window = ContigWindowCls(self.left, self.right)
      self.event_windows.append(new_window)
    else:
      curr_left = self.left
      curr_right = curr_left + window_width - 1
      while (curr_right <= self.right): #{
        new_window = ContigWindowCls(curr_left, curr_right)
        self.event_windows.append(new_window)
        #ExtremeDebugMsg(self, "  Window: %s" % new_window.CoordsString())
        curr_left += 1
        curr_right = curr_left + window_width - 1
      #} end while
    #} end if
    ExtremeDebugMsg(self, "  # Windows: %i" % len(self.event_windows))
    if (0 == len(self.event_windows)):
      LogMsg(self, "Bad region: %i-%i (%ibp), " %
          (self.left, self.right, self.Span()))
      raise R2CGroupError("cannot get any event windows for event region")
    # end if
  #} end def

  def ProcessRead(self, read, unique): #{
    ExtremeDebugMsg(self, "    %s %i-%i Unique: %s" %
      (read.id, read.left, read.right, unique))
    first_window = max(read.left - self.left, 0)
    offset = max(self.right-read.right, 0)+1
    last_window  = len(self.event_windows)-offset
    ExtremeDebugMsg(self, "    First: %i, Last: %i" %
      (first_window, last_window))
    if (first_window <= last_window):
      for window in self.event_windows[first_window:last_window+1]:
        ExtremeDebugMsg(self, "Adding read to window: %i-%i" %
          (window.left, window.right))
        window.AddRead(unique)
      # end for
    # end if
    ExtremeDebugMsg(self, "    [%s]" %
      ",".join(["%2i" % window.support for window in self.event_windows]))
  #} end def

  def InitializeSupport(self, support_string): #{
    self.event_windows = list()
    for window_support in support_string.split(","): #{
      (total, unique) = map(int, window_support.split(":"))
      new_window = ContigWindowCls()
      new_window.support = total
      new_window.support_unique = unique
      self.event_windows.append(new_window)
    #} end for
  #} end def

  def AddSupport(self, support_from): #{
    for (my_window, from_window) in zip(self.event_windows,
        support_from.event_windows): #{
      my_window.AddSupport(from_window)
    #} end for
  #} end def

  def CalculateSupport(self): #{
    self.min_support = None
    self.min_support_unique = None
    total_support = 0
    total_support_unique = 0
    for window in self.event_windows:
      if (None == self.min_support or
          self.min_support > window.support):
        self.min_support = window.support
      # end if
      if (None == self.min_support_unique or
          self.min_support_unique > window.support_unique):
        self.min_support_unique = window.support_unique
      # end if
      total_support += window.support
      total_support_unique += window.support_unique
    # end for
    self.average_support = (
      float(total_support) / float(len(self.event_windows)))
    self.average_support_unique = (
      float(total_support_unique) / float(len(self.event_windows)))
    ExtremeDebugMsg(self, "Support: %i, Unique: %i\nAverage: %.2f, "
      "Unique: %.2f" % (self.min_support, self.min_support_unique,
       self.average_support, self.average_support_unique))
  #} end def

  def HasSupport(self): #{
    for window in self.event_windows: #{
      if (0 < window.support): #{
        return True
      #} end if
    #} end for
    return False
  #} end def

  def SupportString(self): #{
    data_list = [
      self.member_id,
      ",".join(["%i:%i" % (window.support, window.support_unique) for
        window in self.event_windows])
    ]
    return " ".join(data_list)
  #} end def

  def DebugString(self): #{
    data_list = [
      self.member_id,
    ]
    if (None != self.ctg_id): #{
      data_list.append(self.ctg_id)
    #} end if
    if (None != self.left and None != self.right): #{
      data_list.append("%i-%i" % (self.left, self.right))
    #} end if
    if (None != self.event_windows): #{
      data_list.append(",".join(["%i:%i" % (window.support,
        window.support_unique) for window in self.event_windows]))
    #} end if
    return " ".join(data_list)
  #} end def
#} end class

class ContigWindowCls: #{
  def __init__(self, left=None, right=None): #{
    self.left = left
    self.right = right
    self.support = 0
    self.support_unique = 0
  #} end def

  def AddRead(self, unique):
    self.support += 1
    if (unique):
      self.support_unique += 1
    # end if
  # end def

  def AddSupport(self, from_window): #{
    self.support += from_window.support
    self.support_unique += from_window.support_unique
  #} end def

  def CoordsString(self): #{
    return "%i-%i" % (self.left, self.right)
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class R2CGroupError(MyError): #{
  pass
#} end class
