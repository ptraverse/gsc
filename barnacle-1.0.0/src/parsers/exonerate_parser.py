"""
exonerate_parser.py

Created by Readman Chiu
Edited by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

#import standard modules
import re, os
from string import ljust

# import custom modules
import alignment
from utils.general import ReverseComplement

exon = re.compile(r'\s+Query: (?P<query>.+)\n'
                  r'\s+Target: (?P<target>\S+).*\n'
                  r'\s+Model: (?P<model>\S+)\n'
                  r'.+'
                  r'\s+Raw score: (?P<score>\d+)\n'
                  r'\s+Query range: (?P<qstart>\d+) -> (?P<qend>\d+)'
                  r'\s+Target range: (?P<tstart>\d+) -> (?P<tend>\d+)'
                  )
ryo = re.compile(r'ryo (.+)\n')
gff = re.compile(r'--- START OF GFF DUMP ---.+?--- END OF GFF DUMP ---', re.S)
coord = re.compile(r'(\S+):(\d+)-(\d+)')
pairwise = re.compile(r'(?P<pairwise>[^\n]*\d+\s+: .+ :\s+\d+.+\s*\d+\s+: .+ :\s+\d+)', re.S)
vulgar = re.compile(r'vulgar: (.+)\n')

def split(f): #{
    sep = "C4 Alignment:"
    lines = ""
    for line in f: #{
        if line[:len(sep)] == sep: #{
            yield lines
            lines = ""
        lines += line

    yield lines

def parse(file, filters, splice_motif_file=None, noline=True, log_info=None): #{
    "given batch output file, parse results into alignment objects"
    
    #print "parsing exonerate"

    #extracts splice motives
    splice_motives = None

    alignments = []
    
    f = open(file, 'r')

    group = []

    skip_query = None
    for record in split(f): #{
        if record[:len('Command')] == 'Command': #{
            continue
        
        m = exon.search(record)

        #don't create alignment object if regex doesn't work
        if not m or not m.groups or len(m.groups()) < 1: #{
            continue

        align = alignment.Alignment(method="exonerate")
        for (field, value) in m.groupdict().iteritems(): #{
            #make sure the extracted field in an attribute of Alignment
            if field in dir(align): #{
                #print field, value

                if field in ('qstart', 'qend', 'tstart', 'tend', 'score'): #{
                    value = int(value)

                #qstart always different by 1
                if field == 'qstart': #{
                    value += 1

                setattr(align, field, value)

        if skip_query and skip_query == align.query: #{
            continue

        #pairwise alignment
        p = pairwise.search(record)
        if p and p.group('pairwise'): #{
            align.pairwise = p.group('pairwise')
        
        #if ryo output
        ryos = ryo.search(record)
        if ryos: #{
            #print align.query,ryos.group(1)
            for fv in ryos.group(1).split(' '): #{
                (field, value) = fv.split(':')

                if field in dir(align): #{
                    setattr(align, field, value)
        
        if align.tstart > align.tend: #{
            align.target_strand = False
        else:
            align.target_strand = True

        #fix tstart for reverse strand
        if not align.target_strand: #{
            align.tend += 1
        else:
            align.tstart += 1

        #match len = qend - qstart + 1; if ryo reports it, overwrite
        align.match_len = int(align.qend) - int(align.qstart) + 1

        #get blocks
        if align.model.lower() == "est2genome": #{
            align.blocks = get_blocks(record, "exon")
            #have to be done after getting pairwise-alignment
            set_orient(align, splice_motives)
        else:
            align.blocks = get_blocks(record, "similarity")

        #get query blocks
        vulgars = vulgar.search(record)
        if vulgars: #{
            align.query_blocks = get_query_blocks(vulgars.group(1), align.qstart)
        
        if not group or align.query == group[-1].query: #{
            if not group: #{
                align.rank = 1
            elif float(group[-1].score) == float(align.score): #{
                align.rank = group[-1].rank
            else:
                align.rank = int(group[-1].rank) + 1
                
            group.append(align)
        else:
            align.rank = 1
            group = [align]

        #do filtering at beginning to speed up
        if filters: #{
            if 'qlen' in filters: #{
                if not(align.query_len and int(align.query_len) >= int(filters['qlen'])): #{
                    skip_query = align.query
                    continue

            if 'unique' in filters: #{
                #same query i.e same group
                if len(group) > 1 and float(align.score) == float(group[-2].score): #{
                    #first hit would have been already in group, take it out
                    skip_query = align.query
                    if alignments and alignments[-1].query == align.query: #{
                        del alignments[-1]
                    group = []
                    continue

            if 'count' in filters: #{
                if len(group) > int(filters['count']): #{
                    skip_query = align.query
                    continue
            
            if 'bestn' in filters: #{
                if len(group) > 0 and group[-1].rank and int(group[-1].rank) > int(filters['bestn']): #{
                    skip_query = align.query
                    continue
            
            if 'match' in filters: #{
                if align.match_len and align.query_len and float(align.match_len)/float(align.query_len)*100 < float(filters['match']): #{
                    continue

            if 'identity' in filters: #{
                if not (align.identity and float(align.identity) >= float(filters['identity'])): #{
                    continue

            if 'target' in filters: #{
                skip = True
                #just chromosome
                for target in filters['target']: #{
                    if not ':' in target and not '-' in target: #{
                        if align.target.lower() == target.lower(): #{
                            skip = False 

                if skip: #{
                    continue

        skip_query = None

        #keep original keep for output
        if not noline: #{
            align.record = record
            
        alignments.append(align)

    f.close()

    return alignments

def get_blocks(record, method): #{
    "get blocks from gff (vulgar requires more calculations)"

    gffs = gff.findall(record)
    blocks = {}
    if gffs: #{
        #both targetgff and querygff in file, use targetgff(has "Query" in it)
        if len(gffs) > 1: #{
            wanted = [g for g in gffs if "Query" in g][0]
        else:
            wanted = gffs[0]

        if method == "similarity": #{
            a = re.compile("similarity\s+\d+\s+\d+")
            blks = a.findall(wanted)
            blocks = [re.split('\s+', blk)[1:3] for blk in blks]
        elif method == "exon": #{
            e = re.compile("^.+\sexon\s.+$", re.M)
            exons = e.findall(wanted)
            blocks = [re.split('\s+', exon)[3:5] for exon in exons]
            
        return blocks
    else:
        return False

def get_query_blocks(vulgar, qstart): #{
    exons = re.split("I 0 \d+", vulgar)

    p = re.compile('[MG] (\d+)')

    block_sizes = []
    for exon in exons: #{
        qs = p.findall(exon)
        size = 0
        for q in qs: #{
            size += int(q)
        block_sizes.append(size)

    blocks = []

    for size in block_sizes: #{
        if not blocks: #{
            start = int(qstart)
        else:
            start = blocks[-1][1] + 1

        end = start + int(size) - 1

        blocks.append([start, end])

    return blocks

def set_orient(align, splice_motives): #{
    if align.blocks > 1: #{
        (query, aln, target, qstart, qend, tstart, tend, strand) = preprocess_pairwise(align)
        align.splice_sites = get_splice_sites(query, target, align.target_strand)
        align.orient = get_orient(align.splice_sites, align.target_strand, splice_motives)

def get_snvs(align): #{
    sf = alignment.SNVfinder(align)
    snvs = []

    (query, aln, target, qstart, qend, tstart, tend, strand) = preprocess_pairwise(align)
    blocks = get_pairwise_blocks(query, target, aln, qstart, tstart, strand)

    for block in blocks: #{
        snvs.extend(sf.find_snv(block[0],block[1],block[2],block[3],block[4],block[5]))

    return snvs

#pairwise processing
seq_pat = re.compile(r'(?P<start>\d+)\s+: (?P<seq>.+) :\s+(?P<end>\d+)')
intron_pat = re.compile(r'\s{2}[\<\>]+\s+[^\<\>]+\s+[\<\>]+\s{2}')
intron_len_pat = re.compile(r'(\d+)')

def preprocess_pairwise(align): #{
    lines = [l for l in align.pairwise.split("\n") if re.search('\S', l)]

    #concatenate lines together
    query = aln = target = ""
    qstart = qend = tstart = tend = None
    for i in xrange(0, len(lines), 3): #{
        q = seq_pat.search(lines[i])
        query += q.group('seq')

        #right pad alignment line with spaces if necessary
        aln_line = lines[i+1][q.start('seq'):q.end('seq')]
        if len(q.group('seq')) > len(aln_line): #{
            aln_line = ljust(aln_line, len(q.group('seq')))
        aln += aln_line
            
        t = seq_pat.search(lines[i+2])
        target += t.group('seq')

        if i == 0: #{
            qstart = q.group('start')
            tstart = t.group('start')

    #capture end positions when exiting loop
    qend = q.group('end')
    tend = t.group('end')

    if int(tstart) < int(tend): #{
        strand = True
    else:
        strand = False

    return (query, aln, target, qstart, qend, tstart, tend, strand)
 

def join_lines(lines): #{
    seq_pat = re.compile(r'(?P<start>\d+)\s+: (?P<seq>.+) :\s+(?P<end>\d+)')

    query = aln = target = ""
    for i in xrange(0, len(lines), 3): #{
        q = seq_pat.search(lines[i])
        query += q.group('seq')

        #right pad alignment line with spaces if necessary
        aln_line = lines[i+1][q.start('seq'):q.end('seq')]
        if len(q.group('seq')) > len(aln_line): #{
            aln_line = ljust(aln_line, len(q.group('seq')))
        aln += aln_line
            
        t = seq_pat.search(lines[i+2])
        target += t.group('seq')

        if i == 0: #{
            qstart = q.group('start')
            tstart = t.group('start')

    return (query, aln, target)

def get_splice_sites(query, target, strand): #{
    intron_pat = re.compile(r'\s{2}[\<\>]+\s+[^\<\>]+\s+[\<\>]+\s{2}')
    intron_len_pat = re.compile(r'(\d+)')

    splice_sites = []

    if intron_pat.findall(query): #{
        introns = intron_pat.finditer(query)

        for intron in introns: #{
            start = intron.start()
            end = intron.end()
            start_bases = target[start:start+2].lower()
            end_bases = target[end-2:end].lower()
            motif = start_bases + end_bases

            if not strand: #{
                motif = ReverseComplement(motif)
                
            splice_sites.append(motif)

    return splice_sites
                          

def get_orient(splice_sites, strand, splice_motives): #{
    if not splice_motives: #{
        return None
    
    counts = {"forward":0, "backward":0, "unknown":0}

    motives = {'forward':[], 'backward':[]}
    for motif in splice_motives.keys(): #{
        motives['forward'].append(motif.lower())
        motives['backward'].append(ReverseComplement(motif).lower())

    #motives = {'forward':['gtag'], 'backward':['ctac']}

    for ss in splice_sites: #{
        orient = "unknown"

        if ss in motives['forward']: #{
            orient = "forward"
        elif ss in motives['backward']: #{
            orient = "backward"

        counts[orient] += 1

    orient = None
    if counts['forward'] > 0 and counts['backward'] == 0 and counts['forward'] > counts['unknown']: #{
        orient = 'forward'
    elif counts['backward'] > 0 and counts['forward'] == 0 and counts['backward'] > counts['unknown']: #{
        orient = 'backward'

    return orient

def get_pairwise_blocks(query, target, aln, qstart, tstart, strand): #{
    blocks = []

    if intron_pat.findall(query): #{
        introns = intron_pat.finditer(query)
        start_idx = 0
        exon_qstart = int(qstart)
        exon_tstart = int(tstart)

        for intron in introns: #{
            m = intron_len_pat.search(aln[intron.start():intron.end()])
            intron_len = int(m.group())

            exon_query = query[start_idx:intron.start()]
            exon_aln = aln[start_idx:intron.start()]
            exon_target = target[start_idx:intron.start()]
            #adjust end coordinate by subtracting number of deleted bases
            exon_dels = re.findall('-', exon_query)
            target_dels = re.findall('-', exon_target)

            exon_qend = exon_qstart + len(exon_query) - 1 - len(exon_dels)

            if strand: #{
                exon_tend = exon_tstart + len(exon_query) - 1 - len(target_dels)
            else:
                exon_tend = exon_tstart - len(exon_query) + 1 + len(target_dels)

##            print exon_query, exon_qstart, exon_qend, len(exon_query)
##            print exon_aln
##            print exon_target, exon_tstart, exon_tend

            block = [[exon_query, exon_target], exon_aln, [exon_qstart, exon_tstart], True, strand, "exonerate"]
            blocks.append(block)

            start_idx = intron.end()
            exon_qstart = exon_qend + 1
            if strand: #{
                exon_tstart = exon_tend + 1 + intron_len
            else:
                exon_tstart = exon_tend - 1 - intron_len

        exon_query = query[start_idx:len(query)]
        exon_aln = aln[start_idx:len(query)]
        exon_target = target[start_idx:len(query)]
        #adjust end coordinate by subtracting number of deleted bases
        exon_dels = re.findall('-', exon_query)
        target_dels = re.findall('-', exon_target)

        exon_qend = exon_qstart + len(exon_query) - 1 - len(exon_dels)

        if strand: #{
            exon_tend = exon_tstart + len(exon_query) - 1 - len(target_dels)
        else:
            exon_tend = exon_tstart - len(exon_query) + 1 + len(target_dels)

##        print exon_query, exon_qstart, exon_qend, len(exon_query)
##        print exon_aln
##        print exon_target, exon_tstart, exon_tend

        block = [[exon_query, exon_target], exon_aln, [exon_qstart, exon_tstart], True, strand, "exonerate"]
        blocks.append(block)

    else:
        block = [[query, target], aln, [qstart, tstart], True, strand, "exonerate"]
        blocks.append(block)

    return blocks
    

