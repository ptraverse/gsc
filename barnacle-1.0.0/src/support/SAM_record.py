#!/usr/bin/env python
# encoding: utf-8
"""
SAM_record.py

Created by Rod Docking on 2010-10-08.
Edited by Lucas Swanson
Copyright (c) 2010 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

import sys
import os
import string
import re
import subprocess
from utils.general import ReverseComplement

cigar_pattern = re.compile('\d+\D')

class SAM_record:
    """A simple class for holding a single SAM record"""

    def __init__(self, line):

        # Raw line
        self.raw_sam = line
        c = line.split('\t')

        # Check for header lines
        if c[0] in ["@HD", "@SQ", "@RG", "@PG", "@CO"]:
            self.header = True
            self.header_list = c
            return
        else:
            self.header = False

        # Remaining SAM fields
        self.qname = c[0]
        self.flag = int(c[1])
        self.rname = c[2]
        self.pos = int(c[3])
        self.mapq = int(c[4])
        self.cigar = c[5]
        self.mrnm = c[6]
        self.mpos = int(c[7])
        self.isize = int(c[8])
        self.seq = c[9]
        self.qual = c[10]
        self.opt_fields = c[11:]

        # Convenience fields
        self.read_length = None
        if None != re.search(r"^[acgt]+$", c[9].lower()):
            self.read_length = len(c[9])
        self.strand = "Reverse" if (self.flag & 0x10) else "Forward"
        if self.first_read():
            self.read = 1
        elif self.second_read():
            self.read = 2
        else:
            self.read = None

        # Store a few of the more common optional fields (some of these may be BWA-specific)
        self.edit_distance = None
        for field in self.opt_fields:
            if "NM" in field:
                self.edit_distance = int(field.split(":")[-1])

    def __repr__(self):
        return "<SAM_record: {0}>".format(self.qname)

    def __str__(self):
        return """qname (read):\t{0} {1} \
                \nflag:\t\t{2} \
                \nrname (pos):\t{3} {4} \
                \nmapq:\t\t{5} \
                \ncigar:\t\t{6} \
                \nmrnm (pos):\t{7} {8} \
                \nisize:\t\t{9} \
                \nopt_fields:\t{10}\n""".format(
                    self.qname, self.read, self.flag, self.rname, self.pos,
                    self.mapq, self.cigar, self.mrnm, self.mpos, self.isize,
                    self.opt_fields)

    # Flag checks
    def paired_in_sequencing(self):
        return True if (self.flag & 0x1) else False
    def proper_pair(self):
        return True if (self.flag & 0x2) else False
    def q_unmapped(self):
        return True if (self.flag & 0x4) else False
    def m_unmapped(self):
        return True if (self.flag & 0x8) else False
    def reverse_strand(self):
        return True if (self.flag & 0x10) else False
    def mate_reverse_strand(self):
        return True if (self.flag & 0x20) else False
    def first_read(self):
        return True if (self.flag & 0x40) else False
    def second_read(self):
        return True if (self.flag & 0x80) else False
    def not_primary(self):
        return True if (self.flag & 0x100) else False
    def platform_qc_fail(self):
        return True if (self.flag & 0x200) else False
    def pcr_duplicate(self):
        return True if (self.flag & 0x400) else False

    # Utility functions
    def print_fastq(self, re_reverse = True):
        """Print a fastq record from a SAM line"""

        if self.reverse_strand() and re_reverse:
            tmp_seq = ReverseComplement(self.seq)
            tmp_qual = self.qual[::-1]
        else:
            tmp_seq = self.seq
            tmp_qual = self.qual

        return "@{0}/{1}\n{2}\n+\n{3}\n".format(
            self.qname, self.read, tmp_seq, tmp_qual)

    def print_sam(self):
        """Dump the record as SAM"""

        fields = [self.qname, self.flag, self.rname, self.pos, self.mapq,
            self.cigar, self.mrnm, self.mpos, self.isize, self.seq, self.qual]
        fields.extend(self.opt_fields)
        str_fields = [str(c) for c in fields]

        return "{0}".format("\t".join(str_fields))

    def cigar_list(self):
        """Split the cigar field into a list"""

        return re.findall(cigar_pattern, self.cigar)

    def qual_list(self, conversion = 33):
        """Split the quality string into a list of qualities"""

        return [ (ord(qual)-conversion) for qual in self.qual ]

    def perfect_match(self, edit_distance_threshold = 0):
        """ Return true if the record is an exact match to it's reference """

        # Filter out unmapped reads
        if self.q_unmapped():
            return False

        # Confirm that CIGAR is read-length + 'M'
        cigar_pattern = r"^[0-9]+M$"
        if None != self.read_length:
            cigar_pattern = "^%iM$" % self.read_length
        if None == re.search(cigar_pattern, self.cigar):
            return False

        if self.edit_distance > edit_distance_threshold:
            return False

        # If we pass all these checks, return True
        return True

    def perfect_softclipped(self, edit_distance_threshold = 0):
        """ Return true if the record is an exact match to it's reference, with soft-clipping """

        # Filter out unmapped reads
        if self.q_unmapped():
            return False

        # Confirm that CIGAR is read-length + 'M' or xSyMzS where x+y = read-length
        cigar_match = re.search(r"((?P<x>[0-9]+)S)?(?P<y>[0-9]+)M((?P<z>[0-9]+)S)?", self.cigar)
        if None == cigar_match:
            return False
        if None != self.read_length:
          try:
              x = int(cigar_match.group("x"))
          except TypeError:
              x = 0
          y = int(cigar_match.group("y"))
          try:
              z = int(cigar_match.group("z"))
          except TypeError:
              z = 0
          if self.read_length != (x + y + z):
              return False

        if self.edit_distance > edit_distance_threshold:
            return False

        # If we pass all these checks, return True
        return True

    def sanitize_qname(self):
        """ Attempt to remove read numbers from query names """

        orig_qname = self.qname
        if "/" in orig_qname[-2:]:
            self.qname = orig_qname[:-2]
            self.read = int(orig_qname[-1:])


if __name__ == '__main__':
    print "Hey, it's SAM_record.py!"
