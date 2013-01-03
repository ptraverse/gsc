#!/bin/bash

# setup_annotations.sh
#
# Created by Lucas Swanson
# Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.

annots_dir=${0//\/setup_annotations.sh}

# get Ensembl 59 gene annotations and sequences files
echo "Getting ensembl 59 gene annotations file..."
wget ftp://ftp.ensembl.org/pub/release-59/gtf/homo_sapiens/Homo_sapiens.GRCh37.59.gtf.gz --no-host-directories --timestamping --no-verbose --directory-prefix ${annots_dir}/ensembl59
if [[ ! -e ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.gtf || ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.gtf -ot ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.gtf.gz ]]
then
  echo "  decompressing Homo_sapiens.GRCh37.59.gtf.gz"
  gunzip -c ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.gtf.gz > ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.gtf
else
  echo "  Homo_sapiens.GRCh37.59.gtf up to date"
fi
echo "Getting ensembl 59 gene sequences file..."
wget ftp://ftp.ensembl.org/pub/release-59/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.59.cdna.all.fa.gz --no-host-directories --timestamping --no-verbose --directory-prefix ${annots_dir}/ensembl59
if [[ ! -e ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.cdna.all.fa || ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.cdna.all.fa -ot ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.cdna.all.fa.gz ]]
then
  echo "  decompressing Homo_sapiens.GRCh37.59.cdna.all.fa.gz"
  gunzip -c ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.cdna.all.fa.gz > ${annots_dir}/ensembl59/Homo_sapiens.GRCh37.59.cdna.all.fa
else
  echo "  Homo_sapiens.GRCh37.59.cdna.all.fa up to date"
fi
echo

# get Ensembl 65 gene sequences file
echo "Getting ensembl 65 gene sequences file..."
wget ftp://ftp.ensembl.org/pub/release-65/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.65.cdna.all.fa.gz --no-host-directories --timestamping --no-verbose --directory-prefix ${annots_dir}
if [[ ! -e ${annots_dir}/Homo_sapiens.GRCh37.65.cdna.all.fa || ${annots_dir}/Homo_sapiens.GRCh37.65.cdna.all.fa -ot ${annots_dir}/Homo_sapiens.GRCh37.65.cdna.all.fa.gz ]]
then
  echo "  decompressing Homo_sapiens.GRCh37.65.cdna.all.fa.gz"
  gunzip -c ${annots_dir}/Homo_sapiens.GRCh37.65.cdna.all.fa.gz > ${annots_dir}/Homo_sapiens.GRCh37.65.cdna.all.fa
else
  echo "  Homo_sapiens.GRCh37.65.cdna.all.fa up to date"
fi
echo

# get UCSC genes, UCSC ensembl genes, RepeatMasker, and SimpleRepeats coordinates files
echo "Getting UCSC annotation files..."
wget -B ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ -i ${annots_dir}/UCSC_annots_to_get.txt --no-host-directories --timestamping --no-verbose --directory-prefix ${annots_dir}
# create UCSC_genes_ref.txt
if [[ ! -e ${annots_dir}/UCSC_genes_ref.txt || ${annots_dir}/UCSC_genes_ref.txt -ot ${annots_dir}/knownGene.txt.gz || ${annots_dir}/UCSC_genes_ref.txt -ot ${annots_dir}/kgXref.txt.gz ]]
then
  echo "  Creating UCSC_genes_ref.txt"
  gunzip -c ${annots_dir}/knownGene.txt.gz | sort > ${annots_dir}/knownGene_sorted.txt
  gunzip -c ${annots_dir}/kgXref.txt.gz | sort > ${annots_dir}/kgXref_sorted.txt
  join ${annots_dir}/knownGene_sorted.txt ${annots_dir}/kgXref_sorted.txt -t $'\t' > ${annots_dir}/UCSC_genes_ref.txt
  rm ${annots_dir}/knownGene_sorted.txt ${annots_dir}/kgXref_sorted.txt
else
  echo "  UCSC_genes_ref.txt up to date"
fi
# decompress UCSC_genes_hg19.fa
if [[ ! -e ${annots_dir}/UCSC_genes_hg19.fa || ${annots_dir}/UCSC_genes_hg19.fa -ot ${annots_dir}/knownGeneMrna.txt.gz ]]
then
  echo "  decompressing knownGeneMrna.txt.gz to UCSC_genes_hg19.fa"
  gunzip -c ${annots_dir}/knownGeneMrna.txt.gz | awk '{print ">"$1"\n"$2}' > ${annots_dir}/UCSC_genes_hg19.fa
else
  echo "  UCSC_genes_hg19.fa up to date"
fi
# create ensembl65_ref.txt
if [[ ! -e ${annots_dir}/ensembl65_ref.txt || ${annots_dir}/ensembl65_ref.txt -ot ${annots_dir}/ensGene.txt.gz || ${annots_dir}/ensembl65_ref.txt -ot ${annots_dir}/ensemblToGeneName.txt.gz ]]
then
  echo "  Creating ensembl65_ref.txt"
  gunzip -c ${annots_dir}/ensGene.txt.gz | sort -k2 > ${annots_dir}/ensGene_sorted.txt
  gunzip -c ${annots_dir}/ensemblToGeneName.txt.gz | sort > ${annots_dir}/ensemblToGeneName_sorted.txt
  join -1 2 -2 1 ${annots_dir}/ensGene_sorted.txt ${annots_dir}/ensemblToGeneName_sorted.txt -a 1 > ${annots_dir}/ensembl65_ref.txt
  rm ${annots_dir}/ensGene_sorted.txt ${annots_dir}/ensemblToGeneName_sorted.txt
else
  echo "  ensembl65_ref.txt up to date"
fi
# reformat repeat coordinates files
if [[ ! -e ${annots_dir}/hg19_all_rmsk.coords || ${annots_dir}/hg19_all_rmsk.coords -ot ${annots_dir}/rmsk.txt.gz ]]
then
  echo "  Creating hg19_all_rmsk.coords"
  gunzip -c ${annots_dir}/rmsk.txt.gz | awk '{print $6, $7, $8, $11, $12, $13, $2}' > ${annots_dir}/hg19_all_rmsk.coords
else
  echo "  hg19_all_rmsk.coords up to date"
fi
if [[ ! -e ${annots_dir}/hg19_simple_repeats.coords || ${annots_dir}/hg19_simple_repeats.coords -ot ${annots_dir}/simpleRepeat.txt.gz ]]
then
  echo "  Creating hg19_simple_repeats.coords"
  gunzip -c ${annots_dir}/simpleRepeat.txt.gz | awk '{print $2,$3,$4,"TRF_SimpleTandemRepeat_"$17}' > ${annots_dir}/hg19_simple_repeats.coords
else
  echo "  hg19_simple_repeats.coords up to date"
fi
# create the structural RNAs file
if [[ ! -e ${annots_dir}/structural_RNA_repeats.coords || ${annots_dir}/structural_RNA_repeats.coords -ot ${annots_dir}/hg19_all_rmsk.coords ]]
then
  echo "  Creating structural_RNA_repeats.coords"
  awk '{if ($6 ~ /RNA/) print $1, $2, $3, $4, $5, $6}' ${annots_dir}/hg19_all_rmsk.coords | sed "s/^chr//" > ${annots_dir}/structural_RNA_repeats.coords
else
  echo "  structural_RNA_repeats.coords up to date"
fi
echo

# get the UCSC hg19 two-bit human genome sequence file (for event_simulator)
echo "Getting UCSC hg19.2bit genome sequence file..."
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit --no-host-directories --timestamping --no-verbose --directory-prefix ${annots_dir}
echo

echo "setup_annotations.sh COMPLETE"
