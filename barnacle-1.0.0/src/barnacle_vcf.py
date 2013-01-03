"""
barnacle_vcf.py

Created by William Li
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

import sys, os, time, datetime, commands, re
from optparse import OptionParser
sys.path.append("/projects/transabyss/trans-ABySS/v1.2.4/code")

from parsers.candidate_group_parser import CandidateGroupParserCls
from parsers.two_bit import TwoBitFileCls

def get_CHROM(member):
	breakpoint1 = member.breakpointA.ToString()
	breakpoint2 = member.breakpointB.ToString()
	
	breakpoint1 = breakpoint1.partition(':')
	breakpoint2 = breakpoint2.partition(':')
	
	chrom1 = breakpoint1[0]
	chrom1 = chrom1[3:]
	
	chrom2 = breakpoint2[0]
	chrom2 = chrom2[3:]
	
	return chrom1, chrom2

def get_POS(member):
	breakpoint1 = member.breakpointA.ToString()
	breakpoint2 = member.breakpointB.ToString()
	
	overlap = int(member.meta_fields['ctg_overlap'])
	
	position1 = breakpoint1.partition(':')
	position2 = breakpoint2.partition(':')
	
	position1 = position1[2]
	position2 = position2[2]
	
	position1 = position1.partition('(')[0]
	position2 = position2.partition('(')[0]
	
	if overlap > 0:
	    if breakpoint2.endswith('(down)'):
	    	position2 = int(position2)-overlap
	    
	    if breakpoint1.endswith('(down)'):
	    	position1 = int(position1)-overlap  		
	
	return str(position1), str(position2)

def get_REF(member, refseq):
	pos1, pos2 = get_POS(member)
	
	chrom1, chrom2 = get_CHROM(member)
	seq1 = refseq.GetSequence(chrom1, int(pos1), int(pos1))
	seq2 = refseq.GetSequence(chrom2, int(pos2), int(pos2))
	
	return seq1, seq2
	
def get_ALT(member, refseq):
	chrom1, chrom2 = get_CHROM(member)
	
	pos1, pos2 = get_POS(member)
	
	base1, base2 = get_REF(member, refseq)
	
	breakpoint1 = member.breakpointA.ToString()
	breakpoint2 = member.breakpointB.ToString()
	
	overlap = int(member.meta_fields['ctg_overlap'])
	if overlap > 0:
	    if breakpoint1.endswith('(down)') and breakpoint2.endswith('(down)'):
	    	pos1 = str(int(pos1) + overlap)
	    	pos2 = str(int(pos2) + overlap)
	
	#both breakpoints at either the start or the end of the contig region
	if breakpoint1.endswith('(up)') and breakpoint2.endswith('(up)'):
	    if overlap > 0:
	    	pos1 = int(pos1)+overlap
	    	pos2 = int(pos2)+overlap
	    		    	    	    	    	
	    alt1 = '['+chrom2+':'+str(pos2)+'['+base1
	    alt2 = '['+chrom1+':'+str(pos1)+'['+base2
	    	
	    return alt1, alt2
	elif breakpoint1.endswith('(down)') and breakpoint2.endswith('(down)'):	    	    
	    alt1 = base1+']'+chrom2+':'+str(pos2)+']'
	    alt2 = base2+']'+chrom1+':'+str(pos1)+']'
	    
	    return alt1, alt2
	
	#one breakpoint is at the start of the contig region and other breakpoint is at the end
	if breakpoint1.endswith('(up)'):
	    alt1 = ']'+chrom2+':'+pos2+']'+base1
	    alt2 = base2+'['+chrom1+':'+pos1+'['
	else:
	    alt1 = base1+'['+chrom2+':'+pos2+'['
	    alt2 = ']'+chrom1+':'+pos1+']'+base2
	    
	return alt1, alt2
	
def get_QUAL():
	return '.'
	
def get_FILT():
	return 'PASS'
	
def get_INFO(member, id1, id2):
	overlap = str(member.meta_fields['ctg_overlap'])
	
	svtype = 'FND'
	
	dp = int(member.avg_read_to_ctg_unique)
	
	if int(overlap) > 0:
	    return 'SVTYPE='+svtype+';MATEID='+str(id2)+'b;CIPOS=0,'+overlap+';SR='+str(dp)+';CTG=', svtype+';MATEID='+str(id1)+'a;CIPOS=0,'+overlap+';SR='+str(dp)+';CTG='
	else:
	    return 'SVTYPE='+svtype+';MATEID='+str(id2)+'b;SR='+str(dp)+';CTG=', svtype+';MATEID='+str(id1)+'a;SR='+str(dp)+';CTG='

#header output method	
def write_header(GIN_user, GIN_pass, LIMS_user, LIMS_pass, refseq_flag, library, filetype_flag, out_file, contig=None):
	#file format
	out_file.write('##fileformat=VCFv4.1\n')
	    
	#file date
	out_file.write('##filedate='+time.strftime("%Y%m%d")+'\n')
	    
	#tcga version
	out_file.write('##tcgaversion=1.0\n')
	    
	#genome reference; need to use URL
	if refseq_flag == 'hg19':
	    out_file.write('##reference=<ID=hg19,Source=http://www.bcgsc.ca/downloads/genomes/Homo_sapiens/hg19/1000genomes/bwa_ind/genome/\n')
	elif refseq_flag == 'hg18':
	    out_file.write('##reference=<ID=hg18,Source=http://www.bcgsc.ca/downloads/genomes/Homo_sapiens/hg18/bwa_ind/genome/tcga_ref/>\n')
	elif refseq_flag == 'mm9':
	    out_file.write('##reference=<ID=mm9,Source=/projects/transabyss/trans-ABySS/annotations/mm9/201107/genome.fa>\n')
	    
	#contig assembly tags, need to use URL
	out_file.write('##assembly='+contig+'\n')
	    
	#center
	out_file.write('##center="BCGSC"\n')
	    
	#phasing
	out_file.write('##phasing=none\n')
	
	info_format = {
		'svtype':'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n',
		'mateid':'##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakends">\n',
		'event':'##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of breakend event">\n',
		'cipos':'##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">\n',
		'svlen':'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
		'end':'##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n',
		'inv':'##ALT=<ID=INV,Description="Inversion">\n',
		'del':'##ALT=<ID=DEL,Description="Deletion">\n',
		'duptan':'##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">\n',
		'sr':'##INFO=<ID=SR,Number=1,Type=Integer,Description="Spanning reads">\n',
		'dp':'##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n',
		'CTG':'##INFO=<ID=CTG,Number=.,Type=String,Description="Contig ID">\n'
		}
		
	fusion = ['svtype', 'mateid', 'cipos', 'CTG', 'sr']
	duplication = ['svtype', 'duptan', 'CTG', 'sr']
	
	if filetype_flag == 'fusion':
	    for item in fusion:
	        if item in info_format:
	            out_file.write(info_format[item])
	elif filetype_flag == 'itd' or filetype_flag == 'ptd':
	    for item in duplication:
	    	if item in info_format:
	    	    out_file.write(info_format[item])
	
	sample_info = commands.getoutput("python get_tcga_sample_info.py --username "+GIN_user+" --password "+GIN_pass+" --LIMS_user "+LIMS_user+" --LIMS_pass "+LIMS_pass+" --library "+library)    
	
	#sample info
	sample_info = sample_info.split(',')
	patient = sample_info[0]
		    
	sample_id = sample_info[1]
	sample_desc = sample_info[2]
	platform = sample_info[3]
	accession = sample_info[4]
	out_file.write('##SAMPLE=<ID='+sample_id+',Individual='+patient+',Description="'+sample_desc+'",Platform='+platform+',Accession='+accession+'>\n')
	    
	#pedigree
	out_file.write('##PEDIGREE=<Name_0='+sample_id+'>\n')
	
	fields_line = '#CHROM\t'+'POS\t'+'ID\t'+'REF\t'+'ALT\t'+'QUAL\t'+'FILTER\t'+'INFO\n'
	out_file.write(fields_line)

def create_fusion_dict(output_dir, in_file, refseq, sequence_dict):
	fusion_dict = {}
		
	groupParser = CandidateGroupParserCls(in_file)
	
	id1, id2 = 1, 1
	key1, key2 = 1, 2
	
	ctg_dir = output_dir.strip('vcf')
	contigfa = open(ctg_dir+'fa', "w")
	
	for group in groupParser:
	    member = group.members[0]
	    chrom1, chrom2 = get_CHROM(member)
	    	
	    pos1, pos2 = get_POS(member)
	    	
	    ref1, ref2 = get_REF(member, refseq)
	    	
	    alt1, alt2 = get_ALT(member, refseq)
	    	
	    qual = get_QUAL()
	    filt = get_FILT()
	    	
	    info1, info2 = get_INFO(member, id1, id2)
	    	
	    fusion1 = chrom1+'\t'+pos1+'\t'+str(id1)+'a'+'\t'+ref1+'\t'+alt1+'\t'+qual+'\t'+filt+'\t'+info1
	    fusion2 = chrom2+'\t'+pos2+'\t'+str(id2)+'b'+'\t'+ref2+'\t'+alt2+'\t'+qual+'\t'+filt+'\t'+info2
	    
	    counter = 0
	    for m in group.members:
	    	contig = m.contig_info.ToString()
		contig = contig.replace(':', '_')
		contig = contig.partition('(')[0]
		
		if counter > 0:
	    	    fusion1 += ','+contig
	    	    fusion2 += ','+contig
	    	else:
	    	    fusion1 += contig
	    	    fusion2 += contig
	    	
	    	sequence = sequence_dict['>'+contig]    
	    	contigfa.write('>'+contig+'\n'+sequence+'\n')
	    	    
	    	counter += 1
	    
	    fusion_dict[key1] = fusion1+'\n'
	    fusion_dict[key2] = fusion2+'\n'
	    	
	    id1 += 1
	    id2 += 1
	    
	    key1 += 2
	    key2 += 2
	    	
	return fusion_dict

def dup_CHROM(member):
	breakpoint = member.breakpointA.ToString()
	
	breakpoint = breakpoint.partition(':')
	
	chrom = breakpoint[0]
	chrom = chrom[3:]
	
	return chrom

def dup_POS(member):
	breakpoint = member.breakpointA.ToString()
	breakpoint = breakpoint.partition(':')
		
	position = breakpoint[2]
	position = position.partition('(')[0]
		
	return position
	
def dup_REF(member, refseq):
	pos = dup_POS(member)
	
	chrom = dup_CHROM(member)
	seq = refseq.GetSequence(chrom, int(pos), int(pos))
		
	return seq
	
def dup_ALT():
	alt = '<DUP:TANDEM>'
	
	return alt
	
def dup_INFO(member):
	pos = int(dup_POS(member))
	svtype = 'FND'
	length = len(member.event_seq)
	end = pos + length - 1
	
	info = 'SVTYPE='+svtype+';END='+str(end)
	
	return info
	
def ptd_POS(member):
	breakpoint1 = member.breakpointA.ToString()
	breakpoint2 = member.breakpointB.ToString()
	
	if breakpoint1.endswith('(down)'):
	    end = breakpoint2
	    start = breakpoint1.partition(':')[2]
	else:
	    end = breakpoint1
	    start = breakpoint2.partition(':')[2]
	
	start = start.partition('(')[0]
	
	return start
	
def ptd_INFO(member):
	breakpoint1 = member.breakpointA.ToString()
	breakpoint2 = member.breakpointB.ToString()
	
	if breakpoint1.endswith('(down)'):
	    end = breakpoint2
	    start = breakpoint1.partition(':')[2]
	else:
	    end = breakpoint1
	    start = breakpoint2.partition(':')[2]
	
	start = start.partition('(')[0]
	
	endpos = end.partition(':')[2]
	endpos = endpos.partition('(')[0]
	
	svtype = 'FND'
	
	info = 'SVTYPE='+svtype+';END='+endpos
	
	return info

def create_dup_dict(in_file, refseq, filetype_flag):
	dup_dict = {}
	
	groupParser = CandidateGroupParserCls(in_file)
	
	id1 = 1
	
	for group in groupParser:
	    member = group.members[0]
	    
	    chrom = dup_CHROM(member)
	    	
	    qual = get_QUAL()
	    filt = get_FILT()
	    
	    if filetype_flag == 'itd':	
	    	pos = dup_POS(member)
	    else:
	    	pos = ptd_POS(member)
	    	
	    ref = dup_REF(member, refseq)
	    
	    alt = dup_ALT()
	    
	    if filetype_flag == 'itd':
	    	info = dup_INFO(member)
	    else:
	    	info = ptd_INFO(member)
	    	
	    dp = int(member.avg_read_to_ctg_unique)
	    info += ';SR='+str(dp)+';CTG='

	    dup = chrom+'\t'+pos+'\t.\t'+ref+'\t'+alt+'\t'+qual+'\t'+filt+'\t'+info
	    
	    counter = 0
	    for m in group.members:
	    	contig = m.contig_info.ToString()
		contig = contig.replace(':', '_')
		contig = contig.partition('(')[0]
		
		if counter > 0:
	    	    dup += ','+contig
	    	else:
	    	    dup += contig
	    	    
	    	counter += 1
	    	
	    dup_dict[id1] = dup+'\n'
	    	
	    id1 = id1+1
	    	
	return dup_dict

#parse contig(fa) file and put reads into dictionary with contig ID as key
def parse_fa(fa_file):
	sequence = {}
	contig = None
	for line in open(fa_file, 'r'):
     	    if line[0] == '>':
         	contig = line.split()[0]
     	    else:
         	sequence[contig] = line.rstrip('\n')
        
        return sequence

def out_to_VCF(filetype_flag, in_file, gene_flag, output_dir, contig, library, GIN_user, GIN_pass, LIMS_user, LIMS_pass):
	#reference sequence file
	if gene_flag == 'hg18':
	    refseq = TwoBitFileCls('/projects/transabyss/trans-ABySS/annotations/hg18/200909/hg18.2bit')
	elif gene_flag == 'hg19':
	    refseq = TwoBitFileCls('/projects/transabyss/trans-ABySS/annotations/hg19/201110/hg19.2bit')
	elif gene_flag == 'mm9':
	    refseq = TwoBitFileCls('/projects/transabyss/trans-ABySS/annotations/mm9/201107/mm9.2bit')
	
	out_file = open(output_dir, "w")
	
	write_header(GIN_user, GIN_pass, LIMS_user, LIMS_pass, gene_flag, library, filetype_flag, out_file, contig)
	
	if filetype_flag == 'fusion':
	    sequence_dict = parse_fa(contig)
	    dictionary = create_fusion_dict(output_dir, in_file, refseq, sequence_dict)
	elif filetype_flag == 'itd' or filetype_flag == 'ptd':
	    dictionary = create_dup_dict(in_file, refseq, filetype_flag)
	
	for key in dictionary:
	    out_file.write(dictionary[key])
	
	out_file.close()
	
if __name__ == '__main__':
	parser = OptionParser()
	
	parser.add_option("-f",dest="fusion_filename",default=None,help="directory and name of input fusion file to convert into vcf")
	parser.add_option("-i",dest="itd_filename",default=None,help="directory and name of input itd file to convert into vcf")
	parser.add_option("-p",dest="ptd_filename",default=None,help="directory and name of input ptd file to convert into vcf")
	parser.add_option("-d",dest="output_dir",default=None,help="directory and name of output file to write into")
	parser.add_option("--username",dest="username",default=None,help="GIN username")
    	parser.add_option("--password",dest="password",default=None,help="GIN password")
    	parser.add_option("--LIMS_user",dest="lims_user",default=None,help="LIMS username")
    	parser.add_option("--LIMS_pass",dest="lims_pass",default=None,help="LIMS password")
    	parser.add_option("--library",dest="library",default=None,help="Library name")
    	parser.add_option("-c",dest="contigs_file",default=None,help="directory and name of contig assembly file(*.fa)")
    	parser.add_option("-g",dest="gene_flag",default=None,help="Specify which genome to use")
    	
    	(options, args) = parser.parse_args()
    	
    	if options.fusion_filename:
    	    filetype_flag = 'fusion'
    	    in_file = options.fusion_filename
    	elif options.itd_filename:
    	    filetype_flag = 'itd'
    	    in_file = options.itd_filename
    	elif options.ptd_filename:
    	    filetype_flag = 'ptd'
    	    in_file = options.ptd_filename
    	
	output_dir = '/projects/transabyss/workspace/jira/APA-69/barnacle_test/A05012_barnacle_fus.vcf'
	library = 'A05012'
	in_file = '/genesis/scratch/validations/transcriptome/AML/A05012/Assembly/current/barnacle/ver_1.2/8_predicted_events/A05012.barnacle.fus'
	contig_file = '/genesis/scratch/validations/transcriptome/AML/A05012/Assembly/current/merge/A05012-contigs.fa'
	
	out_to_VCF('fusion', in_file, 'hg19', output_dir, contig_file, library, gin_user, gin_pass, lims_user, lims_pass)





