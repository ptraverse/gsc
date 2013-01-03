"""
psl_parser.py

Created by Readman Chiu
Edited by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

#import standard modules
import sys, os, math, re

# import custom modules
import alignment
from utils.general import ReverseComplement
from utils.messages import DebugMsg

fields = {1:"match",
          2:"mismatch",
          3:'repmatch',
          4:'ncount',
          5:'qnuminsert',
          6:'qbaseinsert',
          7:'tnuminsert',
          8:'tbaseinsert',
          9:"query_strand",
          10:"query",
          11:"query_len",
          12:"qstart",
          13:"qend",
          14:"target",
          15:"target_len",
          16:"tstart",
          17:"tend",
          18:"block_count",
          19:"block_sizes",
          20:"qstarts",
          21:"tstarts"
          }

def parse(file, filters=None, noblocks=False, noline=False, minimum=False, splice_motif_file=None, refseq=None, genome=None, log_info=None): #{
    aligns = []

    prev_query = None
    prev_line = None

    group = []
    cols = []

    splice_motifs = None

    chrom_proper = None
    #if genome: #{
    #    chrom_proper = ucsc_chroms(genome)

    for line in open(file, 'r'): #{
        line = line.rstrip()

        #double lines seen in blat output?
        if line == prev_line: #{
            continue

        #check if line begins with number - not headers
        if line and ord(line[0]) >= ord('0') and ord(line[0]) <= ord('9'): #{
            #DebugMsg(log_info, line)
            if not filters: #{
                aligns.append(create_align(line, noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs, chrom_proper=chrom_proper))
                continue

            cols = line.split("\t")
            if (10 > len(cols)): #{
              print "Error parsing line"
              print line
              raise Exception
            #} end if

            if filters and filters.has_key('exclude') and filters['exclude'].has_key(cols[9]): #{
                continue

            #filtering
            if prev_query and not cols[9] == prev_query: #{
                try:
                    indices = screen(group, filters)
                except ValueError, e:
                    print "Error screening group"
                    print line
                # end try

                #only extract blocks for filtered set
                if indices: #{
                    for idx in indices: #{
                        a = create_align(group[idx], noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs, chrom_proper=chrom_proper)
                        aligns.append(a)

                del group[:]

            group.append(line)

            prev_query = cols[9]
            prev_line = line
            del cols[:]

    #last group
    if group and filters: #{
        try:
            indices = screen(group, filters, log_info=log_info)
        except ValueError, e:
            print "Error screening last group"
            print line
        # end try
        if indices: #{
            for idx in indices: #{
                a = create_align(group[idx], noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs, chrom_proper=chrom_proper)
                aligns.append(a)


    return aligns

def create_align(line, noblocks=False, noline=False, minimum=False, refseq=None, splice_motifs=None, chrom_proper=None): #{
    cols = line.split("\t")

    align = alignment.Alignment(method="blat")
    for field in fields: #{
        if fields[field] in ('query', 'target', 'query_strand', 'mismatch'): #{
            setattr(align, fields[field], cols[field-1])

        if not minimum: #{
            if fields[field] in ["qstart", "tstart"]: #{
                setattr(align, fields[field], int(cols[field-1])+1)
            elif fields[field] in ["qend", "tend"]: #{
                setattr(align, fields[field], int(cols[field-1]))
            elif fields[field] in ('qstarts', 'tstarts', 'block_sizes'): #{
                continue
            else:
                setattr(align, fields[field], cols[field-1])

    if not minimum: #{
        try:
          align.score = calc_score(cols[0], cols[1], cols[4], cols[6])
        except ValueError, e:
          print "Error parsing line"
          #print line
          raise e
        # end try
        align.match_len = int(cols[0]) + int(cols[1]) +  int(cols[2])
        align.identity = calc_identity(int(cols[11])+1, cols[12], int(cols[15])+1, cols[16], cols[4], cols[1], align.match_len)

    if not noline: #{
        setattr(align, 'psl_str', line)

    if not noblocks: #{
        align.blocks = get_blocks(cols[20], cols[18])
        align.query_blocks = get_blocks(cols[19], cols[18], strand=cols[8], qsize=cols[10])

    if (refseq or splice_motifs or chrom_proper) and not re.match('^(chr|scaffold)', align.target, re.IGNORECASE): #{
	align.target = 'chr' + align.target


    #splice sites and orientation
    if align.blocks and refseq: #{
        get_splice_sites(align, refseq)
        if splice_motifs: #{
            align.orient = get_orient(align.splice_sites, splice_motifs)

    if chrom_proper and chrom_proper.has_key(align.target): #{
	chrom = align.target
	if chrom_proper.has_key(align.target): #{
	    chrom = chrom_proper[align.target]
	    align.target = chrom

    return align


def get_blocks(starts, block_sizes, strand=None, qsize=None): #{
    blocks = []
    starts = starts.rstrip(',').split(',')
    blk_sizes = block_sizes.rstrip(',').split(',')
    if strand and qsize and strand == '-': #{
        for i in range(len(starts)): #{
            end = int(qsize)-int(starts[i])
            start = end - int(blk_sizes[i]) + 1
            blocks.append([end, start])
    else:
        for i in range(len(blk_sizes)): #{
            block = [int(starts[i])+1, int(starts[i])+int(blk_sizes[i])]
            blocks.append(block)
        #blocks.sort(lambda x,y: x[0]-y[0])

    return blocks

def screen(group, filters, log_info=None): #{
    if not filters: #{
        return range(len(group))

    keep = []

    scores = {}
    ranks = {}

    for idx in range(len(group)): #{
        cols = group[idx].split("\t")
        try:
          scores[idx] = calc_score(cols[0], cols[1], cols[4], cols[6])
        except ValueError, e:
          print "Error parsing group"
          #print group[idx]
          raise e
        # end try

    indices = scores.keys()
    indices.sort(lambda x,y: scores[y]-scores[x])

    #rank
    rank = 1
    for ii in range(len(indices)): #{
        if 'bestn' in filters and rank > int(filters['bestn']): #{
            break

        idx = indices[ii]

        if ii == 0: #{
            ranks[idx] = 1
        elif scores[idx] == scores[indices[ii-1]]: #{
            ranks[idx] = ranks[indices[ii-1]]
        else:
            rank += 1
            ranks[idx] = rank

    #skips entire group
    #unique means best alignment only has 1 hit
    #if 'unique' in filters and len(indices) > 1 and (scores[indices[0]] - scores[indices[1]]) < 0.10 * scores[indices[0]]: #{
    if 'unique' in filters and len(indices) > 1 and scores[indices[0]] == scores[indices[1]]: #{
        return keep

    count = 1
    for idx in indices: #{
        cols = group[idx].split("\t")

        if filters.has_key('bestn') and ranks[idx] > int(filters['bestn']): #{
            break

        if filters.has_key('count') and count > int(filters['count']): #{
            break

        if filters.has_key('qlen') and int(cols[10]) < int(filters['qlen']): #{
            break

        match_len = int(cols[0]) + int(cols[1]) +  int(cols[2])
        if filters.has_key('match'): #{
            if float(match_len)*100/float(cols[10]) < float(filters['match']): #{
                continue

        if filters.has_key('identity'): #{
            #100% identity, just look at number of mismatches
            if float(filters['identity']) == 100.0 and int(cols[1]) > 0: #{
                continue
            else:
                #DebugMsg(log_info, "calc_identity(int(cols[%s])+1, cols[%s], "
                #  "int(cols[%s])+1, cols[%s], cols[%s], cols[%s], "
                #  "match_len)" % (fields[11], fields[12], fields[15],
                #  fields[16], fields[4], fields[1]))
                identity = calc_identity(int(cols[11])+1, cols[12], int(cols[15])+1, cols[16], cols[4], cols[1], match_len)
                if identity < float(filters['identity']): #{
                    continue

        if filters.has_key('ungapped') and int(cols[17]) > 1: #{
            continue

        count += 1
        keep.append(idx)

    return keep

def calc_score(match, mismatch, qnuminsert, tnuminsert): #{
    "ucsc score"
    return int(match) - int(mismatch) - int(qnuminsert) - int(tnuminsert)

def calc_identity(qstart, qend, tstart, tend, qnuminsert, mismatch, match_len): #{
    "ucsc webblat identity"

    q_alisize = int(qend) - int(qstart)
    t_alisize = int(tend) - int(tstart)
    alisize = min(q_alisize, t_alisize)

    if alisize <= 0: #{
        return 0

    sizediff = q_alisize - t_alisize
    if sizediff < 0: #{
        sizediff = 0

    insert_factor = int(qnuminsert)

    millibad = (1000 * (float(mismatch) + insert_factor + round(3 * math.log(1 + sizediff)))) / float(match_len)

    return round(100.0 - millibad * 0.1, 1)

def get_splice_sites(align, refseq): #{
    if refseq.has_key(align.target): #{
        for i in range(len(align.blocks)-1): #{
            if align.blocks[i+1][0] - align.blocks[i][1] < 5: #{
                ss = "NA"
            else:
                intron = refseq[align.target].seq.tostring()[align.blocks[i][1]:align.blocks[i+1][0]-1]
                ss = intron[:2].lower() + intron[-2:].lower()
            align.splice_sites.append(ss)
            #print i, ss
    #print "splice",align.splice_sites, len(align.splice_sites)

def get_orient(splice_sites, splice_motifs): #{
    "copied from exonerate.py"
    if not splice_motifs: #{
        return None

    counts = {"forward":0, "backward":0, "unknown":0}
    motifs = {'forward':[], 'backward':[]}
    for motif in splice_motifs.keys(): #{
        motifs['forward'].append(motif.lower())
        motifs['backward'].append(ReverseComplement(motif).lower())
    for ss in splice_sites: #{
        orient = "unknown"
        if ss in motifs['forward']: #{
            orient = "forward"
        elif ss in motifs['backward']: #{
            orient = "backward"
        counts[orient] += 1
    orient = None
    if counts['forward'] > 0 and counts['backward'] == 0 and counts['forward'] > counts['unknown']: #{
        orient = '+'
    elif counts['backward'] > 0 and counts['forward'] == 0 and counts['backward'] > counts['unknown']: #{
        orient = '-'
    return orient

if __name__ == '__main__': #{
    parse(sys.argv[1])
