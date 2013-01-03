#! /usr/bin/env python
"""
region.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules

# import custom modules

class EventRegionGroupCls: #{
  def __init__(self, log_info=None): #{
    self.log_info = log_info
    self.left     = None
    self.right    = None
    # regions is a dictionary of regions,
    #   keyed by their left and right locations
    self.regions  = dict()
    # region_keys is a dictionary of keys for the regions dictionary,
    #   keyed by event member ID
    self.region_keys = dict()
    # reads is a dictionary of ReadToContigAlignCls objects,
    #   keyed by read id
    self.reads    = dict()
  #} end def

  def AddRegion(self, new_region, member_id): #{
    new_key = "%i-%i" % (new_region.left, new_region.right)
    # add the mapping from member ID to region key
    if (member_id in self.region_keys): #{
      raise R2CSupportError("Adding region for member %s twice!" % member_id)
    else:
      DebugMsg(self, "Mapping member %s to region %s (%i)" %
        (member_id, new_key, new_region.Span()))
      self.region_keys[member_id] = new_key
    #} end if
    # if the region does not yet exist
    if (new_key in self.regions): #{
      DebugMsg(self, "Region already found: %s" % new_key)
    else:
      DebugMsg(self, "Adding region to group: %s" % new_key)
      # add the new region
      self.regions[new_key] = new_region
      # if this is the first region in the group
      if (1 == len(self.regions)): #{
        # initialize the group based on the current region
        DebugMsg(self, "Initializing region group")
        self.left   = new_region.left
        self.right  = new_region.right
      else:
        # update the group span
        DebugMsg(self, "Updating region group")
        self.left  = min(self.left,  new_region.left)
        self.right = max(self.right, new_region.right)
      #} end if
    #} end if
  #} end def

  def CalculateSupport(self): #{
    if (0 < len(self.reads)): #{
      for region in self.regions.itervalues(): #{
        region.CalculateSupport(self.reads)
      #} end for
    #} end if
  #} end def

  def GetRegionForMember(self, member_id): #{
    region_key = self.region_keys[member_id]
    return self.regions[region_key]
  #} end def
#} end class

class EventRegionCls: #{
  def __init__(self, member, min_overlap, log_info=None): #{
    self.log_info = log_info
    self.InitializeFromMember(member)
    self.min_overlap = min_overlap
    self.min_support = 0
    self.min_support_unique = 0
    self.average_support = 0
    self.average_support_unique = 0
  #} end def

  def InitializeFromMember(self, member): #{
    self.ctg_id = member.contig_info.id
    self.ctg_len = member.contig_info.length
    if (member.gap): #{
      #self.left  = member.meta_fields['gap_start']
      #self.right = member.meta_fields['gap_end']
      self.left  = member.align_info_B.ctg_start
      self.right = member.align_info_B.ctg_end
    else:
      self.left  = \
        min(member.align_info_A.ctg_end, member.align_info_B.ctg_start)
      self.right = \
        max(member.align_info_A.ctg_end, member.align_info_B.ctg_start)
    #} end if
  #} end def

  def BufferedLeft(self): #{
    # widen the region to satisfy the minimum overlap requirement
    buffered_left = self.left - self.min_overlap
    # ensure that the event does not extend past the edges of the contig
    return max(buffered_left, 2)
  #} end def

  def BufferedRight(self): #{
    # widen the region to satisfy the minimum overlap requirement
    buffered_right = self.right + self.min_overlap
    # ensure that the event does not extend past the edges of the contig
    return min(buffered_right, self.ctg_len-1)
  #} end def

  def Span(self): #{
    return (self.right - self.left) + 1
  #} end def

  def BufferedSpan(self): #{
    return (self.BufferedRight() - self.BufferedLeft()) + 1
  #} end def

  def CalculateSupport(self, reads): #{
    DebugMsg(self, "Region: %i-%i (%ibp), Buffered: %i-%i (%ibp)" %
      (self.left, self.right, self.Span(),
       self.BufferedLeft(), self.BufferedRight(), self.BufferedSpan()))
    windows = self.GetEventWindows()
    if (0 == len(windows)): #{
      LogMsg(self, "Bad region: %i-%i (%ibp), " %
        (self.left, self.right, self.Span()) + "Buffered: %i-%i (%ibp)" %
        (self.BufferedLeft(), self.BufferedRight(), self.BufferedSpan()))
      raise R2CSupportError \
        ("cannot get any event windows for event region")
    #} end if
    self.CalculateSupportForWindows(reads, windows)
    self.min_support = None
    self.min_support_unique = None
    total_support = 0
    total_support_unique = 0
    # windows = self.OldGetEventWindows(min_read_length)
    for window in windows: #{
      #(window_support, window_support_unique) = \
      #  self.CalculateSupportForWindow(window, reads)
      if (None == self.min_support or
          self.min_support > window.support):
        self.min_support = window.support
      #} end if
      if (None == self.min_support_unique or
          self.min_support_unique > window.support_unique):
        self.min_support_unique = window.support_unique
      #} end if
      total_support += window.support
      total_support_unique += window.support_unique
    #} end for
    self.average_support = float(total_support) / float(len(windows))
    self.average_support_unique = \
      float(total_support_unique) / float(len(windows))
    DebugMsg(self, "Support: %i, Unique: %i\nAverage: %.2f, Unique: %.2f" %
      (self.min_support, self.min_support_unique,
       self.average_support, self.average_support_unique))
  #} end def

  def GetEventWindows(self): #{
    windows = list()
    window_width = (2 * self.min_overlap) + 1
    if (self.BufferedSpan() < window_width): #{
      new_window = ContigWindowCls(self.BufferedLeft(), self.BufferedRight())
      windows.append(new_window)
    else:
      left = self.BufferedLeft()
      right = left + window_width - 1
      while (right <= self.BufferedRight()): #{
        new_window = ContigWindowCls(left, right)
        windows.append(new_window)
        #DebugMsg(self, "  Window: %i-%i" % new_window)
        left += 1
        right = left + window_width - 1
      #} end while
    #} end if
    DebugMsg(self, "  # Windows: %i" % len(windows))
    return windows
  #} end def

  def CalculateSupportForWindows(self, reads, windows): #{
    for read in reads.itervalues(): #{
      DebugMsg(self, "    %s %i-%i Unique: %s" %
        (read.id, read.left, read.right, read.unique))
      first_window = max(read.left - self.BufferedLeft(), 0)
      offset = max(self.BufferedRight()-read.right, 0)+1
      last_window  = len(windows)-offset
      DebugMsg(self, "    First: %i, Last: %i" % (first_window, last_window))
      if (first_window <= last_window): #{
        for window in windows[first_window:last_window+1]: #{
          #DebugMsg(self, "Adding read to window: %i-%i" %
          #  (window.left, window.right))
          window.AddRead(read.unique)
        #} end for
      #} end if
      DebugMsg(self, "    [%s]" %
        ",".join(["%2i" % window.support for window in windows]))
    #} end for
  #} end def

  def CalculateSupportForWindow(self, window, reads): #{
    support = support_unique = 0
    DebugMsg(self, "  Window: %s-%s" % (window.left, window.right))
    for read in reads.itervalues(): #{
      DebugMsg(self, "    %s %i-%i Unique: %s" %
        (read.id, read.left, read.right, read.unique))
      if (read.left  <= window.left and
          read.right >= window.right): #{
        DebugMsg(self, "    counting read")
        support += 1
        if (read.unique): #{
          support_unique += 1
        #} end if
      else:
        DebugMsg(self, "    NOT counting read")
      #} end if
    #} end for
    DebugMsg(self, "  Window support: %i, Unique: %i" %
      (support, support_unique))
    return (support, support_unique)
  #} end def
#} end class
