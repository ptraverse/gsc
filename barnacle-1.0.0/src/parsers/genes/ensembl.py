"""
ensembl.py

Created by Readman Chiu
Edited by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

import transcript
from optparse import OptionParser
import os, re

#for ensGene.txt from UCSC
fields_a = {1:"name", 2:"chrom", 3:"strand", 4:"txStart", 5:"txEnd",
  6:"cdsStart", 7:"cdsEnd", 8:"exonCount", 9:"exonStarts", 10:"exonEnds",
  12:"alias"}

#for ensGene_ref.txt created in-house
fields_b = {0:"name", 2:"chrom", 3:"strand", 4:"txStart", 5:"txEnd",
  6:"cdsStart", 7:"cdsEnd", 8:"exonCount", 9:"exonStarts", 10:"exonEnds",
  16:"alias"}

def set_fields(file=None, line=None):
    sep = name_field = None
    fields = fields_a
    
    #determine which ensGene file it is
    if file:
        for l in open(file, 'r'):
            line = l
            break
    if line:
        if line[:3].lower() != 'ens':
            fields = fields_a
            sep = "\t"
            name_field = 1
        else:
            fields = fields_b
            sep = " "
            name_field = 0

    return sep, name_field, fields

def parse(file):
    txts = []
    sep, name_field, fields = set_fields(file=file)
    
    for line in open(file, 'r'):
        cols = line.rstrip("\n").split(sep)

        if cols[0]:
            txt = transcript.Transcript(cols[name_field])
            
            for i in range(len(cols)):
                if i in fields:
                    if fields[i] == 'chrom' and cols[i][:3] != 'chr':
                        cols[i] = 'chr' + cols[i]
                    
                    if i <= 10 or i == 16 or i == 12:
                        setattr(txt, fields[i], cols[i])

            exonStarts = cols[9].rstrip(',').split(',')
            exonEnds = cols[10].rstrip(',').split(',')
            txt.exons = []
            for e in range(len(exonStarts)):
                #start+1: seems necessary at least for mouse ensembl file
                txt.exons.append([int(exonStarts[e])+1, int(exonEnds[e])])

            #calculate transcript length for coverage
            for exon in txt.exons:
                txt.length += int(exon[1]) - int(exon[0]) + 1
            #print txt.name, txt.exonCount, txt.length, txt.exons[0]
            txts.append(txt)

    return txts

def parse_line(line):
    sep, name_field, fields = set_fields(line=line)

    cols = line.rstrip("\n").split(sep)
    if sep and len(cols) > 1:
        txt = transcript.Transcript(cols[name_field])
        
        for i in range(len(cols)):
            if i in fields:
                if fields[i] == 'chrom' and cols[i][:3] != 'chr':
                        cols[i] = 'chr' + cols[i]
                        
                if i <= 10 or i == 16 or i == 12:
                    setattr(txt, fields[i], cols[i])

        exonStarts = cols[9].rstrip(',').split(',')
        exonEnds = cols[10].rstrip(',').split(',')
        txt.exons = []
        for e in range(len(exonStarts)):
            txt.exons.append([int(exonStarts[e])+1, int(exonEnds[e])])

        #calculate transcript length for coverage
        for exon in txt.exons:
            txt.length += int(exon[1]) - int(exon[0]) + 1

        return txt

    return None

def index(input, output):
    sep, name_field, fields = set_fields(file=input)
    
    indices = {}
    data_file = os.path.abspath(input)
    line_num = 1
    for line in open(input, 'r'):
        cols = line.rstrip().split(sep)

        start = int(int(cols[4])/1000)
        end = int(int(cols[5])/1000)
        target = cols[2]
        
        if not re.match('^(chr|scaffold)', target, re.IGNORECASE):
            target = 'chr' + target
        
        #print cols[0],target,start,end
        for n in range(start,end+1):
            index = ':'.join((target,str(n)))
            value = str(line_num)

            if not indices.has_key(index):
                indices[index] = [value]
            else:
                indices[index].append(value)

        line_num += 1

    index_file = open(output, 'w')
    for index in sorted(indices.keys()):
        index_file.write(' '.join((index, ','.join(indices[index]))) + "\n")

def output(txts, outfile):
    fields = fields_a

    list_size = int(fields.keys()[-1])+1

    field_idx = {}
    for idx, field in fields.iteritems():
        if field in ('exonStarts', 'exonEnds', 'exonCount'):
            field_idx[field] = idx

    out = open(outfile, 'w')
    for i in range(len(txts)):
        txt = txts[i]
        
        data = []
        for idx in range(list_size):
            data.append('NA')
        
        for idx, field in fields.iteritems():
            try:
                value = getattr(txt, field)
            except AttributeError:
                continue
            else:
                data[idx] = str(value)

        data[0] = str(i)

        data[field_idx['exonStarts']] = ','.join([str(int(i[0])-1) for i in txt.exons])
        data[field_idx['exonEnds']] = ','.join([str(i[1]) for i in txt.exons])
        data[field_idx['exonCount']] = str(len(txt.exons))
        
        out.write('\t'.join(data) + '\n')
        
    out.close()
    
if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()

    if options.index:
        index(args[0], options.index)
