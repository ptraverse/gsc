"""
alignment.py

Created by Readman Chiu
Edited by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

#import standard modules
import sys, re

# import custom modules
from utils.general import ReverseComplement

class Alignment: #{

    line = None
    method = None
    query = target = None
    query_len = target_len = None
    #strand = Bolean
    target_strand = None
    query_strand = None
    qstart = qend = tstart = tend = None
    blocks = None
    splice_sites = None
    query_blocks = None
    match_len = matches = mismatches = None
    identity = None
    score = None
    rank = None
    model = None
    pairwise = None
    psl_str = None
    #transcript orientation inferred from Exonerate alignment
    orient = None
    contig = None

    def __init__(self, method): #{
        self.method = method
        self.blocks = self.query_blocks = self.splice_sites = []

    def details(self): #{
        if self.blocks: #{
            return "%s %s-%s %s %s-%s %s" % (self.query, self.qstart, self.qend, self.target, self.tstart, self.tend, len(self.blocks))
        else:
            return "%s %s-%s %s %s-%s" % (self.query, self.qstart, self.qend, self.target, self.tstart, self.tend)
        
    def gff(self, type): #{
        if not self.blocks: #{
            print "no blocks", self.query, self.target
            return ""
        
        #chromosome
        if 'chr' in self.target: #{
            chrom = self.target
        else:
            chrom = 'chr' + self.target

        #contig
        contig = self.query.split(" ")[0]
        if self.contig and self.contig.num: #{
            contig = str(self.contig.num)

        #strand
        if self.target_strand: #{
            strand = "+"
        else:
            strand = "-"

        #gene orientation
        orient = "?"
        if self.orient: #{
            if self.orient == 'forward': #{
                orient = "+"
            elif self.orient == 'backward': #{
                orient = "-"

        #coverage
        cov = None
        if self.contig: #{
            cov = str(self.contig.normalized_coverage())

        if cov: #{
            group = ":".join([self.target, str(self.tstart), "%s%s" % (strand, orient), contig, str(cov)])
        else:
            group = ":".join([self.target, str(self.tstart), "%s%s" % (strand, orient), contig])

        #do this after recording strand info in group
        if orient != "?": #{
            strand = orient

        out = ""
        
        for i in range(len(self.blocks)): #{
            qstart, qend = self.query_blocks[i]
            tstart, tend = self.blocks[i]
            
            splice_site = ""
            if len(self.blocks) > 1 and i < len(self.blocks)-1 and len(self.blocks): #{
                splice_site = self.splice_sites[i]

            #hide query block and splice sites in "program" field
            annot = ":".join([str(qstart), str(qend), splice_site])

            out += "\t".join([chrom, annot, type, tstart, tend, ".", strand, ".", group]) + "\n"

        return out

            
    def psl(self): #{
        "output ucsc psl format"

        if self.psl: #{
            psl_str = self.psl_str

            cols = self.psl_str.split("\t")
            #make sure target has chr in it for UCSC
            if not re.match('^(chr|scaffold)', cols[13], re.IGNORECASE): #{
                cols[13] = 'chr' + cols[13]
            
            if self.contig: #{
                cov = self.contig.normalized_coverage()
                if cov: #{
                    contig = ":".join([str(self.contig.num), str(cov)])
                    cols[9] = contig

            #gene orientation
            if self.orient: #{
                if self.orient == 'forward': #{
                    orient = "+"
                elif self.orient == 'backward': #{
                    orient = "-"
                cols[8] = orient
                
            psl_str = "\t".join(cols)
            
            return psl_str + "\n"

    def exon(self): #{
        "output original exonerate record"

        if self.record: #{
            return self.record + "\n"
        
        
    def correct_blocks(self, splice_motifs, target_seq, query_seq):                  
        if not self.splice_sites or not splice_motifs: #{
            return False
	
	#fix fake insertion
	self.correct_unaligned(query_seq, target_seq)
	
	#fix single gaps
	self.correct_single_gaps(splice_motifs, target_seq)
	
	#fix neighboring gaps
	self.correct_neighbor_gaps(splice_motifs, target_seq)
	
	
    def correct_single_gaps(self, splice_motifs, target_seq): #{
	gaps = {}
	
	for i in range(len(self.blocks)-1): #{
            ss = self.splice_sites[i]

            if ss and not splice_motifs.has_key(ss) and not splice_motifs.has_key(ReverseComplement(ss).lower()): #{
                if abs(self.query_blocks[i+1][0] - self.query_blocks[i][1]) == 1: #{
                    gaps[i] = 0
		    
        if gaps: #{
	    target_blocks = self.blocks[:]
	    query_blocks = self.query_blocks[:]
	    splice_sites = self.splice_sites[:]
	    
            gap_indices = gaps.keys()
            gap_indices.sort(lambda x,y: x-y)
	    
            for i in gap_indices: #{
                tblock1 = target_blocks[i][:]
                tblock2 = target_blocks[i+1][:]
                qblock1 = query_blocks[i][:]
                qblock2 = query_blocks[i+1][:]
                splice_site = self.fix_single_gap(tblock1, tblock2, qblock1, qblock2, splice_motifs, target_seq, self.query_strand)

                if splice_site: #{
                    sys.stderr.write("%s changed blocks %s %s to %s %s\n" % (self.query, self.target, target_blocks[i], tblock1, splice_site))
                    sys.stderr.write("%s changed blocks %s %s to %s %s\n" % (self.query, self.target, target_blocks[i+1], tblock2, splice_site))
                    target_blocks[i] = tblock1
                    target_blocks[i+1] = tblock2
                    query_blocks[i] = qblock1
                    query_blocks[i+1] = qblock2
                    splice_sites[i] = splice_site

            if target_blocks != self.blocks: #{
                self.blocks = target_blocks[:]
                self.query_blocks = query_blocks[:]
                self.splice_sites = splice_sites
                if not self.mismatch or int(self.mismatch) == 0: #{
                    self.mismatch = 1
		    	
    def correct_neighbor_gaps(self, splice_motifs, target_seq): #{
	gaps = {}
        
        for i in range(len(self.blocks)-1): #{
            ss = self.splice_sites[i]

            if ss and not splice_motifs.has_key(ss) and not splice_motifs.has_key(ReverseComplement(ss).lower()): #{
                if abs(self.query_blocks[i+1][0] - self.query_blocks[i][1]) == 1: #{
                    gaps[i] = 0
                    
        if gaps: #{
	    target_blocks = self.blocks[:]
	    query_blocks = self.query_blocks[:]
	    splice_sites = self.splice_sites[:]
	    
            gap_indices = gaps.keys()
            gap_indices.sort(lambda x,y: x-y)
            
            #fix by moving exon and then shuffle
            replaced = {}
            replaced_ordered = []
            for i in range(len(gap_indices)-1): #{
                i1 = gap_indices[i]
                i2 = i1 + 1

                if gaps.has_key(i2) and gaps[i2] != 'fixed' and i2 + 1 < len(target_blocks): #{
                    tblock1 = target_blocks[i1][:]
                    tblock2 = target_blocks[i2][:]
                    tblock3 = target_blocks[i2+1][:]
                    qblock1 = query_blocks[i1]
                    qblock2 = query_blocks[i2][:]
                    qblock3 = query_blocks[i2+1][:]

                    splice_site = self.fix_neighbor_gaps(tblock1, tblock2, tblock3, qblock1, qblock2, qblock3, splice_motifs, target_seq, self.query_strand)

		    if splice_site: #{
                        idx = ' '.join((str(i1), str(i2), str(i2+1)))
                        replaced[idx] = tblock1, tblock3, qblock1, qblock3, splice_site
                        gaps[i1] += 1
                        gaps[i2] += 1
                        replaced_ordered.append(idx)

            #make sure delete from back to front
            replaced_ordered.reverse()
            for indices in replaced_ordered: #{
                new_blocks = replaced[indices]

                ok = True
                for index in indices.split(' '): #{
                    if gaps.has_key(int(index)) and gaps[int(index)] > 1: #{
                        ok = False
                        break

                if ok: #{
                    idx = [int(i) for i in indices.split(' ')]
                    sys.stderr.write("%s changed blocks %s %s to %s %s\n" % (self.query, self.target, self.blocks[idx[0]], new_blocks[0], new_blocks[-1]))
                    sys.stderr.write("%s changed blocks %s %s to %s %s\n" % (self.query, self.target, self.blocks[idx[2]], new_blocks[1], new_blocks[-1]))
                    sys.stderr.write("%s removed block %s %s\n" % (self.query, self.target, self.blocks[idx[1]]))
                    target_blocks[idx[0]] = new_blocks[0]
                    target_blocks[idx[2]] = new_blocks[1]
                    query_blocks[idx[0]] = new_blocks[2]
                    query_blocks[idx[2]] = new_blocks[3]
                    splice_sites[idx[0]] = new_blocks[-1]
                    del target_blocks[idx[1]]
                    del query_blocks[idx[1]]
                    del splice_sites[idx[1]]

            if target_blocks != self.blocks: #{
                self.blocks = target_blocks[:]
                self.query_blocks = query_blocks[:]
                self.splice_sites = splice_sites

                if not self.mismatch or int(self.mismatch) == 0: #{
                    self.mismatch = 1
	

    def correct_unaligned(self, query_seq, target_seq, max_diff=5): #{
	"""tries to fix unaligned query sequence because of snvs"""
	
	#print self.query_blocks
	#print self.blocks
	#print self.splice_sites
	
	expand= {}
	for i in range(len(self.blocks)-1): #{
	    if abs(self.query_blocks[i+1][0] - self.query_blocks[i][1]) != 1 and abs(self.blocks[i+1][0] - self.blocks[i][1]) != 1: #{
		if self.query_strand == '+': #{
		    qseq = query_seq[self.query_blocks[i][1]:self.query_blocks[i+1][0]-1]
		else:
		    qseq = ReverseComplement(query_seq[self.query_blocks[i+1][0]:self.query_blocks[i][1]-1]).lower()
		
		tseq = target_seq[self.blocks[i][1]:self.blocks[i+1][0]-1].lower()
		
		gap_seq = tseq[len(qseq):]
		ss_start = gap_seq[:2] + gap_seq[-2:]
		gap_seq = tseq[:-1*len(qseq)]
		ss_end = gap_seq[:2] + gap_seq[-2:]
			
		#if gap in target sequence greater than gap in query sequence, method doesn't apply
		#also don't apply if target gap is less than 20bp(arbitrary) of query seq
		#print "test", i, self.blocks[i], self.blocks[i][1], self.blocks[i+1][0]-1, qseq, tseq
		if len(qseq) > len(tseq) or len(tseq) - len(qseq) < 20: #{
		    continue
		
		#find differences
		start_diff = 0
		for j in range(len(qseq)): #{
		    if qseq[j].lower() != tseq[j].lower(): #{
			start_diff += 1
			
		end_diff = 0
		for j in range(len(qseq)): #{
		    if qseq[j].lower() != tseq[-1*len(qseq):][j].lower(): #{
			end_diff += 1
		
		if start_diff < end_diff and start_diff < max_diff: #{
		    expand[i] = ['left', len(qseq)]
		    self.mismatch = start_diff
		elif end_diff < start_diff and end_diff < max_diff: #{
		    self.mismatch = end_diff
		    expand[i] = ['right', len(qseq)]
	
	for (idx, where) in expand.iteritems(): #{
	    if where[0] == 'left': #{
		sys.stderr.write("%s changed blocks %s\n" % (self.query, self.blocks[idx]))
		if self.query_strand == '+': #{
		    self.query_blocks[idx][1] += where[1]
		else:
		    self.query_blocks[idx][1] -= where[1]
		self.blocks[idx][1] += where[1]
		self.splice_sites[idx] = ss_start
	    else:
		sys.stderr.write("%s changed blocks %s\n" % (self.query, self.blocks[idx+1]))
		if self.query_strand == '+': #{
		    self.query_blocks[idx+1][0] -= where[1]
		else:
		    self.query_blocks[idx+1][0] += where[1]
		self.blocks[idx+1][0] -= where[1]
		self.splice_sites[idx] = ss_end
		
	#print self.query_blocks
	#print self.blocks
	#print self.splice_sites
		
	
    def fix_single_gap(self, tblock1, tblock2, qblock1, qblock2, splice_motifs, target_seq, query_strand): #{
        min_size = 10
        max_size = 10000
        shuffle_sizes = [-2, -1, 1, 2]

        tgap = [tblock1[1]+1, tblock2[0]-1]
        tsize = tgap[1] - tgap[0] + 1

        if tsize < min_size or tsize > max_size: #{
            return False
        
        artefact = False
        ambiguous = False
        splice_site = None
        for size in shuffle_sizes: #{
	    #if shuffling size is greater than or equal to block size, can't shuffle
	    if abs(size) >= (tblock1[1] - tblock1[0]) or abs(size) >= (tblock2[1] - tblock2[0]): #{
		continue
	    
            coord = tgap[0] + size, tgap[1] + size
            gap_seq = target_seq[coord[0]-1:coord[1]]
            
            ss = gap_seq[:2] + gap_seq[-2:]

            if splice_motifs.has_key(ss.lower()) or splice_motifs.has_key(ReverseComplement(ss).lower()): #{
                if not artefact: #{
                    artefact = True
                else:
                    ambiguous = True

                if not ambiguous: #{
                    tblock1[1] += size
                    tblock2[0] += size

                    if query_strand == '+': #{
                        qblock1[1] += size
                        qblock2[0] += size
                    else:
                        qblock1[1] -= size
                        qblock2[0] -= size

                    splice_site = ss

        if artefact and not ambiguous: #{
            return splice_site
        else:
            return False

    def fix_neighbor_gaps(self, tblock1, tblock2, tblock3, qblock1, qblock2, qblock3, splice_motifs, target_seq, query_strand): #{
        tgap1 = [tblock1[1]+1, tblock2[0]-1]
        tgap2 = [tblock2[1]+1, tblock3[0]-1]

        tgap1_size = tgap1[1] - tgap1[0] + 1
        tgap2_size = tgap2[1] - tgap2[0] + 1

        tblock2_size = tblock2[1] - tblock2[0] + 1

        max_shuffle_size = 10
        if tblock2_size > max_shuffle_size: #{
            return False
        
        small_cutoff = 10
        big_cutoff = 20

        big = None
        
        if tgap1_size < small_cutoff and tgap2_size > big_cutoff: #{
            big = 'right'

        elif tgap2_size < small_cutoff and tgap1_size > big_cutoff: #{
            big = 'left'

        if not big: #{
            return False

        if big == 'right': #{
            tblock3[0] -= tblock2_size
            qblock3[0] = qblock2[0]
        else:
            tblock1[1] += tblock2_size
            qblock1[1] = qblock2[1]

        tgap = [tblock1[1]+1, tblock3[0]-1]
        gap_seq = target_seq[tgap[0]-1:tgap[1]]
        ss = gap_seq[:2] + gap_seq[-2:]

        merged_gaps = None
        splice_site = None

        if splice_motifs.has_key(ss.lower()) or splice_motifs.has_key(ReverseComplement(ss).lower()): #{
            splice_site = ss
        else:
            #print "testing11", tblock1, tblock3, qblock1, qblock3
            splice_site = self.fix_single_gap(tblock1, tblock3, qblock1, qblock3, splice_motifs, target_seq, query_strand)
            #print "testing22", tblock1, tblock3, qblock1, qblock3

        return splice_site
