#!/bin/bash

gfagz=$1
refprefix=$2
hap1=$3
hap2=$4
threads=$5

./vcf_extract_phased.sh $gfagz $refprefix $threads

gfa=$(basename $gfagz .gz)
og=$(basename $gfa .gfa).og
fasta=$(basename $gfa .gfa).fa

odgi build -g $gfa -o $og -t $threads
odgi paths -i $og -f >$fasta
samtools faidx $fasta
samtools faidx $fasta $(odgi paths -i $og -L | grep "^$hap1") >$fasta.$hap1.fa
samtools faidx $fasta $(odgi paths -i $og -L | grep "^$hap2") >$fasta.$hap2.fa

#pirs simulate 


