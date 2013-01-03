#! /usr/bin/env python
"""
contig_with_alignments.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from operator import attrgetter, itemgetter

# import custom modules
from utils.error import MyError
from utils.general import SetupMainClass, AddChr
from utils.messages import LogMsg, DebugMsg, ExtremeDebugMsg
from alignment_functions import (CalcAlignOverlap, CalcOverlap, FixAlign,
  ShortAlignString)
from gap_filter import GapFilterCls

class ContigWithAlignmentsCls: #{
  def __init__(self, contig_align_index, aligns, ctg_seq_file,
               gap_out_path, options, log_info):
    SetupMainClass(self, options, log_info=log_info)
    self.align_index = contig_align_index
    self.id = aligns[contig_align_index].query
    self.length = int(aligns[contig_align_index].query_len)
    self.num_aligns_to_contig = 0
    self.gapped_event_found   = "N"
    self.num_gapped_aligns    = 0
    self.single_align_found   = False
    self.perfect_align_found  = False
    self.multi_mapped         = False
    # Create an empty list to hold the best alignments
    # for the current contig
    self.best_aligns = list()
    # Create an empty list to hold the alignment groups
    # for the current contig
    self.align_groups = list()
    self.multi_grouped_aligns = dict()
    # setup the gap filter info
    self.gap_filter = GapFilterCls(ctg_seq_file, gap_out_path,
      self.options, self.log_info)
    self.gapped_psl_lines = list()
    self.aligns_for_gap = None
    self.curr_fract = 0.0
  #} end def

  def SelectAlignments(self, aligns): #{
    #LogMsg(self, "Selecting alignments for contig (trunk)...")
    # iterate over the alignments to the current contig
    for curr_align_index in range(self.align_index, len(aligns)+1): #{
      # stop looking when we are onto the next contig
      if (curr_align_index == len(aligns) or
          aligns[curr_align_index].query != self.id):
        #LogMsg(self, "Moving on to next contig") # DEBUG
        break
      #} end if
      #LogMsg(self, "Curr Align Ctg: %s" %
      #  aligns[curr_align_index].query) # DEBUG

      # if a perfect alignment has been found, we do not need to do anything
      if (self.perfect_align_found): #{
        continue
      #} end if

      # fix the alignment
      curr_align = FixAlign(aligns[curr_align_index])
      #LogMsg(self, "%i) %s %i" %
        #(curr_align_index, curr_align.query, curr_align.score)) # DEBUG
      curr_align.id = curr_align_index
      ExtremeDebugMsg(self, "\n".join(["-------------------------",
        "Curr Align Index: %i/%i" % (curr_align_index, len(aligns)-1),
        "  Align S:%i E:%i T:%s SC:%i PID:%.2f" %
        (curr_align.qstart, curr_align.qend, AddChr(curr_align.target),
         curr_align.score, curr_align.identity)]))

      # if the the gap filter should be used
      if (self.options.check_gap): #{
        # update the group of alignments for the gap filter
        self.UpdateAlignsForGap(curr_align)
      #} end if

      # if not checking for split alignments, or a single alignment has been
      # found for the current contig do not consider other alignments for the
      # contig any further
      if (not self.options.check_split or self.single_align_found): #{
        #LogMsg(self, "Single alignment found") # DEBUG
        continue
      #} end if

      # check whether the current alignment is "perfect" or at least long
      # enough to represent the whole contig on its own
      # (do not use the match_len member of the alignment class)
      (align_len, align_fraction) = GetAlignLengthAndFraction(curr_align)
      if (align_len == curr_align.query_len and not curr_align.qgap and
          curr_align.match == curr_align.query_len and
          100.0 == curr_align.identity): #{
        self.perfect_align_found = True
        DebugMsg(self, "Perfect align: %s" % curr_align.query)
      #} end if
      if ("blat" == curr_align.method.lower()): #{
        match_fraction = float(curr_align.match) / float(curr_align.query_len)
      else:
        match_fraction = float(curr_align.score) / float(align_len)
      #} end if
      #DebugMsg(self, "\n".join(["Checking contig representation...",
      #  "Fraction represented: %.2f" % align_fraction,
      #  "Match fraction: %.2f" % match_fraction]))
      # make the match fraction a little less stringent
      if (align_fraction >= self.options.longest_single_align and
          match_fraction >= (self.options.longest_single_align - 0.05)):
        if (not self.perfect_align_found): #{
          DebugMsg(self, "Single align: %s" % curr_align.query)
          ExtremeDebugMsg(self, "  %s" % ShortAlignString(curr_align))
        #} end if
        #DebugMsg(self, "Fraction represented: %.2f\n" % align_fraction +
        #  "Match fraction: %.2f" % match_fraction)
        # release any best alignments found yet
        del self.best_aligns[:]
        # release any alignment groups found yet
        del self.align_groups[:]
        # set flag, rather than break, so that the rest of the
        #   alignments can be checked with the gap filter
        #LogMsg(self, "Skipping alignments for %s" % self.id) # DEBUG
        self.single_align_found = True
        continue
      #else:
      #  DebugMsg(self, "Not a Single Alignment")
      #} end if

      if (self.options.use_quick_chooser): #{
        #DebugMsg(self, "Using quick chooser...")
        self.UpdateBestAligns(curr_align)
      else:
        #DebugMsg(self, "Using smart chooser...")
        self.UpdateAlignGroups(curr_align)
      #} end if
    #} end for

    # check that we are really done with the current contig
    if (curr_align_index < len(aligns) and
        aligns[curr_align_index].query == self.id):
      msg = "Error grouping alignments for contig %s" % self.id
      raise ContigAlignsError(msg)
    #} end if

    # count the number of alignments to the current contig
    self.num_aligns_to_contig = curr_align_index - self.align_index
    ExtremeDebugMsg(self, "Ctg: %s, Num Aligns: %d" %
      (self.id, self.num_aligns_to_contig))
  #} end def

  def UpdateAlignsForGap(self, curr_align): #{
    if (None == self.aligns_for_gap): #{
      ExtremeDebugMsg(self, "Creating group for gaps")
      self.aligns_for_gap = AlignGroupCls(curr_align,
        self.options, self.log_info)
      # only use score and pid for paring alignments for gap check group
      self.aligns_for_gap.paring['prefer_spliced'] = False
      self.aligns_for_gap.paring['prefer_exons']   = False
      # be slightly looser with group for gap alignments
      #self.aligns_for_gap.paring['max_pid_diff'] += 1.0
      # use gap-specific pid paring parameter
      self.aligns_for_gap.paring['mm_max_pid_diff'] = (
        self.options.mm_max_pid_diff_gap)
    else:
      ExtremeDebugMsg(self, "Adding align to group for gaps")
      self.aligns_for_gap.AddAlign(curr_align)
    #} end if
  #} end def

  def UpdateBestAligns(self, curr_align): #{
    # if no alignments have been selected yet
    if (0 == len(self.best_aligns)): #{
      # initialize the best alignments list with the current alignment
      self.best_aligns.append(curr_align)
      ExtremeDebugMsg(self, "Initializing best alignments...\n%s" %
        ShortAlignString(curr_align))
    # if the best alignments list is not empty,
    # compare the current alignment to the best alignments
    else:
      self.CompareToBestAligns(curr_align)
    #} end if
    #LogMsg(self, "Num Best: %i" % len(best_aligns)) # DEBUG
    #LogMsg(self, "-----") # DEBUG
  #} end def

  def CompareToBestAligns(self, align): #{
    align.add_to_bests = True
    align.in_bests = False
    ExtremeDebugMsg(self, "Curr: %s" % ShortAlignString(align))
    for best_index, best_align in enumerate(self.best_aligns): #{
      ExtremeDebugMsg(self, "Best: %s" % ShortAlignString(best_align))
      # get the overlap fraction between the current alignment and
      # the current best alignment
      (overlap, overlap_fraction) = CalcAlignOverlap(align, best_align)
      ExtremeDebugMsg(self, "Overlap: %i Fraction: %f" %
        (overlap, overlap_fraction))
      # if the overlap fraction is great enough
      if (self.options.min_merge_overlap <= overlap_fraction): #{
        self.CompareAlignToBest(align, best, best_index)
        # if the alignment was not good enough to be processed or used to
        # mark multi-mapping do not compare it to any other best alignments
        if (not (align.process or align.mark_mm)): #{
          # do not compare it to any other best alignments
          break
        #} end if
      #} end if
    #} end for
    if (align.add_to_bests): #{
      ExtremeDebugMsg(self, "Adding alignment to: %s" % AddChr(align.target))
      self.best_aligns.append(align)
    #} end if
  #} end def

  def CompareAlignToBest(self, align, best, best_index): #{
    # compare current alignment to the current best alignment
    self.CompareAligns(align, best_align)
    # if the alignment should be ignored, do not add it to the list
    if (align.process): #{
      if (align.replace_best): #{
        align.add_to_bests = False
        # if the alignment has already replaced a best alignment
        if (align.in_bests): #{
          # remove the current best alignment
          del self.best_aligns[best_index]
        else:
          # mark that the alignment has replaced a best alignment
          align.in_bests = True
          # replace the current best alignment with the current alignment
          self.best_aligns[best_index] = align
        #} end if
        if (align.mark_mm): #{
          ExtremeDebugMsg(self, "Marking new best alignment as multi-mapping")
          align.multi_mapped = True
        #} end if
      else:
        # mark both the current alignment and the current best alignment
        # as multi_mapped
        align.multi_mapped = True
        best_align.multi_mapped = True
      #} end if
    else:
      if (align.mark_mm): #{
        ExtremeDebugMsg(self, "Marking best alignment as multi-mapping")
        best_align.multi_mapped = True
      #} end if
      align.add_to_bests = False
    #} end if
  #} end def

  def CompareAligns(self, align, best): #{
    # if the alignment score is not good enough to use for multi-mapping
    if (self.MMScoreTooLow(align.score, best.score)): #{
      ExtremeDebugMsg(self, "Ignoring alignment. Score: %i Max: %i " %
        (align.score, best.score) + "Fract: %.2f Min: %.2f" %
        (self.curr_fract, self.options.mm_min_score_fract))
      align.process = False
      align.replace_best = False
      align.mark_mm = False
    # if the alignment score is not good enough to keep, ignore it
    elif (self.ScoreTooLow(align.score, best.score)): #{
      ExtremeDebugMsg(self, "Ignoring alignment (marking multi-mapping). "
        "Score: %i Max: %i " % (align.score, best.score) +
        "Fract: %.2f Min: %.2f" % (self.curr_fract,
        self.options.min_score_fract))
      align.process = False
      align.replace_best = False
      align.mark_mm = True
    # if the old best alignment score is not good enough
    # to use for multi-mapping
    elif (self.MMScoreTooLow(best.score, align.score)): #{
      ExtremeDebugMsg(self, "Replacing best with alignment. "
        "Score: %i Max: %i " % (align.score, best.score) +
        "Fract: %.2f Min: %.2f" % (self.curr_fract,
          self.options.mm_min_score_fract))
      align.process = True
      align.replace_best = True
      align.mark_mm = False
    # if the current alignment score is much better than the current best,
    # replace the current best alignment
    elif (self.ScoreTooLow(best.score, align.score)): #{
      ExtremeDebugMsg(self, "Replacing best with alignment "
        "(marking multi-mapping). Score: %i Max: %i " % (align.score,
         best.score) + "Fract: %.2f Min: %.2f" % (self.curr_fract,
        self.options.min_score_fract))
      align.process = True
      align.replace_best = True
      align.mark_mm = True
    # if the current alignment score is good enough to process, but not
    # good enough to replace the current best
    else:
      align.process = True
      align.replace_best = False
      align.mark_mm = True
    #} end if
  #} end def

  def UpdateAlignGroups(self, curr_align): #{
    # if no groups have been created yet
    if (0 == len(self.align_groups)): #{
      # initialize the list with a new alignment group
      ExtremeDebugMsg(self, "Creating initial alignment group...")
      self.align_groups.append(AlignGroupCls(curr_align,
        self.options, self.log_info))
      ExtremeDebugMsg(self, "  Initial group S:%i E:%i Span:%i" %
        (self.align_groups[0].ctg_start, self.align_groups[0].ctg_end,
        self.align_groups[0].Span()))
      #} end if
    else:
      self.CompareToAlignGroups(curr_align)
    #} end if
  #} end def

  def CompareToAlignGroups(self, align): #{
    # create a list to hold the indices of the groups the alignment
    # could belong to
    group_overlaps = list()
    # iterate through the current groups to see if the alignment
    # belongs in any
    #LogMsg(self, "Curr    S:%i E:%i Span:%i" % # DEBUG
        #(align.qstart, align.qend, (align.qend - align.qstart) + 1)) # DEBUG
    for group_num, group in enumerate(self.align_groups): #{
      #LogMsg(self, "Group %i S:%i E:%i Span:%i" % # DEBUG
              #(group_num, group.ctg_start, group.ctg_end, # DEBUG
               #(group.ctg_end - group.ctg_start) + 1)) # DEBUG
      (overlap, overlap_fraction) = CalcGroupOverlap(align, group)
      #LogMsg(self, "  O:%i OF:%f" %
      #  (overlap, overlap_fraction)) # DEBUG
      # if the overlap fraction is great enough
      if (self.options.min_merge_overlap <= overlap_fraction): #{
        # add the group to the list
        group_overlap = (group_num, overlap_fraction)
        group_overlaps.append(group_overlap)
      #} end if
    #} end for
    # if no groups were found
    if (0 == len(group_overlaps)): #{
      # add a new group for the alignment
      ExtremeDebugMsg(self, "Creating new alignment group")
      new_group = AlignGroupCls(align, self.options, self.log_info)
      self.align_groups.append(new_group)
      ExtremeDebugMsg(self, "  New alignment group: S:%i E:%i" %
        (new_group.ctg_start, new_group.ctg_end))
    else:
      # if more than one group was found
      if (1 < len(group_overlaps)): #{
        # sort the groups by overlap value
        group_overlaps.sort(key=itemgetter(1), reverse=True)
        # find the first group with less than the maximum overlap
        index = 0
        while (index < len(group_overlaps) and
               group_overlaps[0] <= group_overlaps[index]):
          index += 1
        #} end while
        # slice off the low-overlap groups
        if (index < len(group_overlaps)): #{
          group_overlaps = group_overlaps[:index]
        #} end if
      #} end if
      # if more than one group remains
      if (1 < len(group_overlaps)): #{
        group_list = [ group_overlap[0] for group_overlap in group_overlaps ]
        self.multi_grouped_aligns[align.id] = group_list
        #msg = ("WARNING: Alignment being added to multiple groups: %s" %
        #  self.id)
        #LogMsg(self, msg) # DEBUG
      #} end if
      for group_num, overlap in group_overlaps: #{
        ExtremeDebugMsg(self, "Adding alignment to group %i S:%i E:%i" %
            (group_num,
             self.align_groups[group_num].ctg_start,
             self.align_groups[group_num].ctg_end))
        self.align_groups[group_num].AddAlign(align)
      #} end for
    #} end if
  #} end def

  def PareAlignmentGroups(self): #{
    for group in self.align_groups: #{
      group.Pare()
      # check whether the contig is "multi-mapped"
      if (group.MultiMapped()): #{
        self.multi_mapped = True
        for align in group.best_aligns: #{
          align.multi_mapped = True
        #} end for
      #} end if
    #} end for
    # REMINDER: pare alignments by exons, and whether they are on the same
    # gene disambiguate any remaining multi-grouped alignments
    if (0 < len(self.multi_grouped_aligns)): #{
      self.DisambiguateMultiGroupedAlignments()
    #} end if
  #} end def

  def DisambiguateMultiGroupedAlignments(self): #{
    ExtremeDebugMsg(self, "-"*40 +
      "\nDisambiguating multi-grouped alignments...")
    for align_id, group_list in self.multi_grouped_aligns.iteritems(): #{
      ExtremeDebugMsg(self, "Disambiguating alignment %s from multiple "
        "groups" % align_id)
      group_count = 0
      for group_id in group_list: #{
        group = self.align_groups[group_id]
        align_found = False
        for index, align in enumerate(group.best_aligns): #{
          if (align_id == align.id): #{
            ExtremeDebugMsg(self, "Multi-grouped alignment to %s found in "
              "group %i" % (AddChr(align.target), group_id))
            align_found = True
            group_count += 1
            break
          #} end if
        #} end for
        if (align_found): #{
          if (1 < group_count): #{
            ExtremeDebugMsg(self, "Removing multi-grouped alignment from "
              "group %i" % group_id)
            #DebugMsg(self, "Before: %i" % len(group.best_aligns))
            del group.best_aligns[index:index+1]
            #DebugMsg(self, "After: %i" % len(group.best_aligns))
          #} end if
        elif (self.log_info['debug']): # DEBUG
          ExtremeDebugMsg(self, "Alignment pared from group %i" % group_id)
        #} end if
      #} end for
      if (self.log_info['debug']): #{
        if (0 == group_count): #{
          ExtremeDebugMsg(self, "Alignment pared from all groups") # DEBUG
          pass
        elif (1 < group_count): #{
          LogMsg(self, "WARNING: alignment still multi-grouped after "
            "paring, arbitrarily removed from all groups but the first")
        #} end if
      #} end if
    #} end for
  #} end def

  # REMINDER: add DisambiguateMultiGroupedAlignment method for contents
  #           of outer loop above
  # REMINDER: add CheckGroupForAlignment method for contents of inner
  #           loop above

  def CheckGappedAlignments(self): #{
    ExtremeDebugMsg(self, "-"*40 + "\nChecking gapped alignments...")
    # reset the gap filter
    self.gap_filter.ResetFilter()
    if (not self.options.use_quick_chooser): #{
      # pare down the alignment group so that
      # only the best alignments remain
      self.aligns_for_gap.Pare()
    #} end if
    num_aligns = len(self.aligns_for_gap.best_aligns)
    ExtremeDebugMsg(self, "  %i alignments to check (+%i extra)" %
      (num_aligns, len(self.aligns_for_gap.mm_aligns)))
    num_aligns += len(self.aligns_for_gap.mm_aligns)
    # if too many alignments were found, do not check any of them
    if (self.options.gap_max_num_aligns < num_aligns): #{
      # get the sequence anyways
      self.gap_filter.GetSequenceForContig(self.id)
      ExtremeDebugMsg(self, "NO GAP CHECK: Too many alignments!")
      return
    #} end if
    # check that the contig is not too long
    if (self.options.gap_max_len < self.length): #{
      LogMsg(self, "NO GAP CHECK: too long contig: %i" %
        self.length)
      return False
    #} end if
    # check whether the contig should be considered multi-mapped
    if (1 < num_aligns): #{
      self.gap_filter.multi_mapped = True
    #} end if
    # iterate through the best alignments for the contig
    for align in self.aligns_for_gap.best_aligns: #{
      # if any gapped alignment events are found
      if(self.gap_filter.FindGappedAlignment(align,
          ("" == self.gap_filter.ctg_seq_id))):
        # record their psl lines to be output
        self.gapped_psl_lines.append(align.psl())
      #} end if
    #} end for
    # increment the number of gapped alignments found
    self.gapped_event_found = self.gap_filter.gapped_event_found
    self.num_gapped_aligns += self.gap_filter.num_gapped_aligns
    ExtremeDebugMsg(self, "# Gaps Found (Contig): %i" % self.num_gapped_aligns)
  #} end def

  def MMScoreTooLow(self, test, compare_to): #{
    self.curr_fract = float(test) / float(compare_to)
    if (self.curr_fract < self.options.mm_min_score_fract): #{
      return True
    #} end if
    return False
  #} end def

  def ScoreTooLow(self, test, compare_to): #{
    self.curr_fract = float(test) / float(compare_to)
    if (self.curr_fract < self.options.min_score_fract): #{
      return True
    #} end if
    return False
  #} end def
#} end class

class AlignGroupCls: #{
  # MEMBERS
  #   aligns:    list of alignments, sorted by score
  #   ctg_start: starting contig coordinate
  #   ctg_end:   ending contig coordinate
  #   span():    length of the contig represented

  def __init__(self, align, options, log_info): #{
    self.best_aligns = [align]
    # mm_aligns is a list of alignments that are good enough to be used for
    # determining whether this group should be marked as multi-mapping, but
    # not good enough to be processed
    self.mm_aligns = list()
    self.ctg_start = min(align.qstart, align.qend)
    self.ctg_end   = max(align.qstart, align.qend)
    # maintain values for paring
    self.paring = {}
    self.paring['maintain']           = options.maintain_pared_groups
    self.paring['mm_min_score_fract'] = options.mm_min_score_fract
    self.paring['mm_max_pid_diff']    = options.mm_max_pid_diff
    self.paring['min_score_fract']    = options.min_score_fract
    self.paring['max_pid_diff']       = options.max_pid_diff
    self.paring['prefer_spliced']     = options.prefer_spliced
    self.paring['prefer_exons']       = options.prefer_exons
    self.max_score     = align.score
    self.max_pid       = align.identity
    self.spliced       = align.spliced
    self.exons         = (0 < len(align.overlap))
    self.log_info      = log_info
  #} end def

  def Span(self): #{
    return (self.ctg_end - self.ctg_start) + 1
  #} end def

  def AddAlign(self, align): #{
    align.add_to_aligns = True
    align.pare_others   = False
    if (self.paring['maintain']): #{
      self.CheckNewAlign(align)
    #} end if
    if (align.add_to_aligns): #{
      if (align.pare_others): #{
        self.max_score = align.score
        self.max_pid   = align.identity
        self.spliced   = align.spliced
        self.exons     = (0 < len(align.overlap))
        old_mm_aligns = self.best_aligns + self.mm_aligns
        self.best_aligns = list()
        self.mm_aligns = list()
        for old_align in old_mm_aligns: #{
          if (self.CheckMMAlign(old_align)):
            self.mm_aligns.append(old_align)
          #} end if
        #} end for
      #} end if
      self.best_aligns.append(align)
    elif (self.CheckMMAlign(align)):
      self.mm_aligns.append(align)
    #} end if
    align_start = min(align.qstart, align.qend)
    align_end   = max(align.qstart, align.qend)
    self.ctg_start = min(self.ctg_start, align_start)
    self.ctg_end   = max(self.ctg_end,   align_end)
  #} end def

  def CheckNewAlign(self, align): #{
    # if the alignment score is not good enough to keep, do not add it
    if (self.ScoreTooLow(align.score, self.max_score)): #{
      ExtremeDebugMsg(self, "Not really adding alignment. Score: %i Max: %i " %
        (align.score, self.max_score) + "Fract: %.2f Min: %.2f" %
        (self.curr_fract, self.paring['min_score_fract']))
      align.add_to_aligns = False
      align.pare_others = False
    # if the alignment score is much better than the current max,
    # replace the rest of the group
    elif (self.ScoreTooLow(self.max_score, align.score)): #{
      ExtremeDebugMsg(self, "Replacing group with alignment. "
        "Score: %i Max: %i " % (align.score, self.max_score) +
        "Fract: %.2f Min: %.2f" % (self.curr_fract,
        self.paring['min_score_fract']))
      align.add_to_aligns = True
      align.pare_others = True
    # if the alignment identity is not good enough to keep, do not add it
    elif (self.PIDTooLow(align.identity, self.max_pid)): #{
      ExtremeDebugMsg(self, "Not really adding alignment. PID: %i Max: %i " %
        (align.identity, self.max_pid) + "Diff: %.2f MaxDiff: %.2f" %
        (self.curr_diff, self.paring['max_pid_diff']))
      align.add_to_aligns = False
      align.pare_others = False
    # if the alignment identity is much better than the current max,
    # replace the rest of the group
    elif (self.PIDTooLow(self.max_pid, align.identity)): #{
      ExtremeDebugMsg(self, "Replacing group with alignment. "
        "PID: %i Max: %i " % (align.identity, self.max_pid) +
        "Diff: %.2f MaxDiff: %.2f" % (self.curr_diff,
        self.paring['max_pid_diff']))
      align.add_to_aligns = True
      align.pare_others = True
    # if the prefer splicing flag is set and the group contains spliced
    # alignments and the alignment is not spliced, do not add it
    elif (self.paring['prefer_spliced'] and
        self.spliced and not align.spliced): #{
      ExtremeDebugMsg(self, "Not really adding alignment. Unspliced")
      align.add_to_aligns = False
      align.pare_others = False
    # if the prefer splicing flag is set and the group does not contain
    # spliced alignments and the alignment is spliced, replace the rest
    # of the group
    elif (self.paring['prefer_spliced'] and
        not self.spliced and align.spliced): #{
      ExtremeDebugMsg(self, "Replacing group with alignment. Spliced")
      align.add_to_aligns = True
      align.pare_others = True
    # if the prefer exons flag is set and the group contains exon-overlapping
    # alignments and the alignment does not overlap any exons, do not add it
    elif (self.paring['prefer_exons'] and
        self.exons and (0 == len(align.overlap))):
      ExtremeDebugMsg(self, "Not really adding alignment. No exons")
      align.add_to_aligns = False
      align.pare_others = False
    # if the prefer exons flag is set and the group does not contain
    # any alignments that overlap exons and the alignment overlaps exons,
    # replace the rest of the group
    elif (self.paring['prefer_exons'] and
        not self.exons and (0 < len(align.overlap))):
      ExtremeDebugMsg(self, "Replacing group with alignment. Exons")
      align.add_to_aligns = True
      align.pare_others = True
    else:
      ExtremeDebugMsg(self, "Updating max values")
      self.max_score = max(self.max_score, align.score)
      self.max_pid = max(self.max_pid, align.identity)
      if (align.spliced): #{
        self.spliced = True
      #} end if
      if (0 < len(align.overlap)): #{
        self.exons = True
      #} end if
      align.add_to_aligns = True
      align.pare_others = False
    #} end if
  #} end def

  def CheckMMAlign(self, align): #{
    # check whether the alignment score is good enough to use this alignment
    # for determining whether this group multi-maps
    if (self.MMScoreTooLow(align.score, self.max_score)): #{
      ExtremeDebugMsg(self, "Discarding MM alignment. Score: %i Max: %i " %
        (align.score, self.max_score) + "Fract: %.2f Min: %.2f" %
        (self.curr_fract, self.paring['mm_min_score_fract']))
      return False
    # check whether the alignment identity is good enough to use this alignment
    # for determining whether this group multi-maps
    elif (self.MMPIDTooLow(align.identity, self.max_pid)): #{
      ExtremeDebugMsg(self, "Discarding MM alignment. PID: %i Max: %i " %
        (align.identity, self.max_pid) + "Diff: %.2f MaxDiff: %.2f" %
        (self.curr_diff, self.paring['mm_max_pid_diff']))
      return False
    #} end if
    return True
  #} end def

  def Pare(self): #{
    if (self.paring['maintain']): #{
      return
    #} end if
    ExtremeDebugMsg(self, "Paring alignment group") # DEBUG
    #LogMsg(self, "---------") # DEBUG
    #for i, align in enumerate(self.best_aligns): # DEBUG
      #LogMsg(self, "%i) Score: %s Spliced:%s PID:%.1f T:%s" %
        #(i, align.score, align.spliced, align.identity, align.target)) # DEBUG
    #} end for
    # pare alignments by score
    self.PareByScore()
    if (1 < len(self.best_aligns)): #{
      # pare alignments by percent identity
      self.PareByPID()
    #} end if
    if (self.paring['prefer_spliced'] and 1 < len(self.best_aligns)): #{
      # pare alignments by whether they are spliced
      self.PareBySplicing()
    #} end if
    if (self.paring['prefer_exons'] and 1 < len(self.best_aligns)): #{
      # pare alignments by whether they overlap exons
      self.PareByExons()
    #} end if
  #} end def

  def PareByScore(self): #{
    ExtremeDebugMsg(self, "BY SCORE") # DEBUG
    # sort alignments by score
    self.best_aligns.sort(key=attrgetter('score'), reverse=True)
    # so the maximum score will now be the score of the first alignment
    self.max_score = self.best_aligns[0].score
    # find the first alignment with a score too low for processing
    first_mm_index = 0
    while (first_mm_index < len(self.best_aligns) and
           not self.ScoreTooLow(self.best_aligns[first_mm_index].score,
             self.max_score)):
      first_mm_index += 1
    #} end while
    # find the first alignment with a score too low for multi-mapping
    first_pare_index = first_mm_index
    while (first_pare_index < len(self.best_aligns) and
           not self.MMScoreTooLow(self.best_aligns[first_pare_index].score,
             self.max_score)):
      first_pare_index += 1
    #} end while
    # slice off the low-scoring alignments
    self.mm_aligns = self.best_aligns[first_mm_index:first_pare_index]
    if (first_mm_index < len(self.best_aligns)): #{
      self.best_aligns = self.best_aligns[:first_mm_index]
    #} end if
    for i, align in enumerate(self.best_aligns): # DEBUG
      ExtremeDebugMsg(self,
        "%i) Score: %s Spliced:%s Exons:%i PID:%.1f T:%s" % # DEBUG
        (i, align.score, align.spliced, len(align.overlap),
         align.identity, AddChr(align.target))) # DEBUG
    #} end for
    for i, align in enumerate(self.mm_aligns): # DEBUG
      ExtremeDebugMsg(self,
        "mm%i) Score: %s Spliced:%s Exons:%i PID:%.1f T:%s" % # DEBUG
        (i, align.score, align.spliced, len(align.overlap),
         align.identity, AddChr(align.target))) # DEBUG
    #} end for
  #} end def

  def PareByPID(self): #{
    ExtremeDebugMsg(self, "BY PID") # DEBUG
    # sort alignments by percent identity
    self.best_aligns.sort(key=attrgetter('identity'), reverse=True)
    # so the maximum identity will now be the identity of the first alignment
    self.max_pid = self.best_aligns[0].identity
    # find the first alignment with a percent identity too low for processing
    first_mm_index = 0
    while (first_mm_index < len(self.best_aligns) and
           not self.PIDTooLow(self.best_aligns[first_mm_index].identity,
             self.max_pid)): #{
      first_mm_index += 1
    #} end while
    # find the first alignment with a percent identity too low
    # for multi-mapping
    old_mm_aligns = self.best_aligns[first_mm_index:]
    if (0 < len(self.mm_aligns)): #{
      old_mm_aligns.extend(self.mm_aligns)
      old_mm_aligns.sort(key=attrgetter('identity'), reverse=True)
    #} end if
    first_pare_index = 0
    while (first_pare_index < len(old_mm_aligns) and
           not self.MMPIDTooLow(old_mm_aligns[first_pare_index].identity,
             self.max_pid)): #{
      first_pare_index += 1
    #} end while
    # slice off the low-identity alignments
    self.mm_aligns = old_mm_aligns[:first_pare_index]
    if (first_mm_index < len(self.best_aligns)): #{
      self.best_aligns = self.best_aligns[:first_mm_index]
    #} end if
    #for i, align in enumerate(self.best_aligns): # DEBUG
      #LogMsg(self,
      #  "%i) Score: %s Spliced:%s Exons:%i PID:%.1f T:%s" % # DEBUG
      #  (i, align.score, align.spliced, len(align.overlap),
      #   align.identity, align.target)) # DEBUG
    #} end for
  #} end def

  def PareBySplicing(self): #{
    ExtremeDebugMsg(self, "BY SPLICING") # DEBUG
    # create an empty list to hold the splice alignments
    spliced_aligns = list()
    new_mm_aligns = list()
    # get the spliced alignments
    for align in self.best_aligns: #{
      if (align.spliced): #{
        spliced_aligns.append(align)
      else:
        new_mm_aligns.append(align)
      #} end if
    #} end for
    # if any spliced alignments were found, discard the unspliced alignments
    if (0 < len(spliced_aligns)): #{
      self.best_aligns = spliced_aligns
      self.mm_aligns.extend(new_mm_aligns)
    #} end if
    #for i, align in enumerate(self.best_aligns): # DEBUG
      #LogMsg(self,
      #  "%i) Score: %s Spliced:%s Exons:%i PID:%.1f T:%s" % # DEBUG
      #  (i, align.score, align.spliced, len(align.overlap),
      #   align.identity, align.target)) # DEBUG
    #} end for
  #} end def

  def PareByExons(self): #{
    ExtremeDebugMsg(self, "BY EXONS") # DEBUG
    # create an empty list to hold the exon-overlapping alignments
    exon_aligns = list()
    new_mm_aligns = list()
    # get the exon-overlapping alignments
    for align in self.best_aligns: #{
      if (0 < len(align.overlap)): #{
        exon_aligns.append(align)
      else:
        new_mm_aligns.append(align)
      #} end if
    #} end for
    # if any exon-overlapping alignments were found,
    # discard the alignments that do not overlap exons
    if (0 < len(exon_aligns)): #{
      self.best_aligns = exon_aligns
      self.mm_aligns.extend(new_mm_aligns)
    #} end if
    for i, align in enumerate(self.best_aligns): # DEBUG
      ExtremeDebugMsg(self,
        "%i) Score: %s Spliced:%s Exons:%i PID:%.1f T:%s" % # DEBUG
        (i, align.score, align.spliced, len(align.overlap),
         align.identity, AddChr(align.target))) # DEBUG
    #} end for
  #} end def

  def ScoreTooLow(self, test, compare_to): #{
    self.curr_fract = float(test) / float(compare_to)
    if (self.curr_fract < self.paring['min_score_fract']): #{
      return True
    #} end if
    return False
  #} end def

  def MMScoreTooLow(self, test, compare_to): #{
    self.curr_fract = float(test) / float(compare_to)
    if (self.curr_fract < self.paring['mm_min_score_fract']): #{
      return True
    #} end if
    return False
  #} end def

  def MMPIDTooLow(self, test, compare_to): #{
    self.curr_diff = compare_to - test
    if (self.curr_diff >= self.paring['mm_max_pid_diff']): #{
      return True
    #} end if
    return False
  #} end def

  def PIDTooLow(self, test, compare_to): #{
    self.curr_diff = compare_to - test
    if (self.curr_diff >= self.paring['max_pid_diff']): #{
      return True
    #} end if
    return False
  #} end def

  def MultiMapped(self): #{
    if (1 < (len(self.best_aligns) + len(self.mm_aligns))): #{
      return True
    #} end if
    return False
  #} end def
#} end class

def GetAlignLengthAndFraction(align): #{
  align_len = (align.qend - align.qstart) + 1
  align_fract = float(align_len) / float(align.query_len)
  return (align_len, align_fract)
#} end def

def CalcGroupOverlap(align, group): #{
  return CalcOverlap(align.qstart, align.qend, group.ctg_start, group.ctg_end)
#} end def

#### EXCEPTION CLASSES ####
class ContigAlignsError(MyError): #{
  """Exception raised for errors encountered by the ContigWithAlignments class"""
  pass
#} end class
