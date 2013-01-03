"""
refGene.py

Created by Readman Chiu
Edited by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

import transcript
from optparse import OptionParser
import os, subprocess

fields = {1:"name", 2:"chrom", 3:"strand", 4:"txStart", 5:"txEnd",
  6:"cdsStart", 7:"cdsEnd", 8:"exonCount", 9:"exonStarts", 10:"exonEnds",
  12:"alias"}

def parse(file):
    txts = []
    
    ff = open(file, 'r')
    
    for line in ff.readlines():
        cols = line.rstrip("\n").split("\t")

        if cols[0]:
            txt = transcript.Transcript(cols[1])
            
            for i in range(len(cols)):
                if i in fields:
                    if i < 9 or i == 12:
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
                
            txts.append(txt)
        
    ff.close()
    
    return txts

def parse_line(line):
    cols = line.rstrip("\n").split("\t")

    if cols[0]:
        txt = transcript.Transcript(cols[1])
            
        for i in range(len(cols)):
            if i in fields:
                if i < 9 or i == 12:
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
    indices = {}
    data_file = os.path.abspath(input)
    line_num = 1
    for line in open(input, 'r'):
        cols = line.rstrip().split("\t")

        start = int(int(cols[4])/1000)
        end = int(int(cols[5])/1000)
        target = cols[2]

        if not 'chr' in target:
            target = 'chr' + target

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

def get_peptide(txt_id, pep_file):
    process = subprocess.Popen(["grep", txt_id, pep_file], shell=False, stdout=subprocess.PIPE)
    lines = process.communicate()[0].rstrip("\n").split('\n')
    
    for line in lines:
        cols = line.split('\t')
        if cols[0] == txt_id:
            print cols[1]
        return cols[1]
        
if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()

    if options.index:
        index(args[0], options.index)
