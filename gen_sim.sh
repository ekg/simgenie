#!/bin/bash

gfagz=$1
refprefix=$2
sample=$3
hap1=$4
hap2=$5
threads=$6

echo "extracting phased VCF from the GFA"

../vcf_extract_phased.sh $gfagz $refprefix $threads

prefix=$(basename $gfagz .gfa.gz)
vcf=$prefix.vcf
toplevelvcf=$prefix.toplevel.vcf
truthvcf=$prefix.$sample.truth.vcf

## XXX this is hard-coded to assume a format specific to how we extract the phase above

vcffilter -f 'LV = 0' $vcf | sed s/$refprefix\_$refprefix#/$refprefix#/ >$toplevelvcf.1
vcfkeepsamples $toplevelvcf.1 $(vcfsamplenames $vcf | grep -v "^$sample") | vcffixup - | vcffilter -f 'AC > 0' >$toplevelvcf
vcfkeepsamples $toplevelvcf.1 $sample | vcffixup - | vcffilter -f 'AC > 0' >$truthvcf
rm -f $toplevelvcf.1

gfa=$prefix.gfa
og=$prefix.og
fasta=$prefix.fa

echo "extract reference and selected haplotypes into FASTA"

odgi build -g $gfa -o $og -t $threads
odgi paths -i $og -f >$fasta
samtools faidx $fasta
samtools faidx $fasta $(odgi paths -i $og -L | grep "^$refprefix")  >$fasta.ref.fa
samtools faidx $fasta.ref.fa
samtools faidx $fasta $(odgi paths -i $og -L | grep "^$hap1") >$fasta.$hap1.fa
samtools faidx $fasta $(odgi paths -i $og -L | grep "^$hap2") >$fasta.$hap2.fa
#cat $fasta.$hap1.fa $fasta.$hap2.fa >$fasta.$hap1+$hap2.fa

echo "simulate reads (with pirs) from the diploid FASTA"

pirs simulate -d -l 100 -x 30 -t $threads -d -o $prefix.pirs \
     -B /home/erik/pirs/Profiles/Base-Calling_Profiles/humNew.PE100.matrix.gz \
     -I /home/erik/pirs/Profiles/InDel_Profiles/phixv2.InDel.matrix \
     -G /home/erik/pirs/Profiles/GC-depth_Profiles/humNew.gcdep_150.dat \
     $fasta.$hap1.fa $fasta.$hap2.fa

simfq=$prefix.sim.fq

#cat $prefix.pirs/*.fq | pigz -p $threads >$simfq
cat $prefix.pirs/*.fq >$simfq

echo "apply pangenie to infer the underlying diplotype (estimated_n = $estimated_n)"
#estimated_n=$(odgi paths -i $og -L | cut -f -2 -d'#' | sort | uniq | wc -l)  # this seems to have marginal effects
#time PanGenie -i $simfq -r $fasta.ref.fa -v $toplevelvcf -k 31 -x $estimated_n -o $gfa.pangenie -e 30000000
time PanGenie -i $simfq -r $fasta.ref.fa -v $toplevelvcf -k 31 -o $gfa.pangenie -e 30000000

echo "use rtg eval to compare the truth and called VCFs"

sdf=$fasta.ref.sdf
rm -rf $sdf
rtg format -o $sdf $fasta.ref.fa
bgzip -f $truthvcf && tabix -p vcf $truthvcf.gz
bgzip -f $gfa.pangenie_genotyping.vcf && tabix -p vcf $gfa.pangenie_genotyping.vcf.gz
bgzip -f $gfa.pangenie_phasing.vcf && tabix -p vcf $gfa.pangenie_phasing.vcf.gz
rm -rf $gfa.pangenie_genotyping.vcf.eval
rm -rf $gfa.pangenie_phasing.vcf.eval
rtg vcfeval -b $truthvcf.gz -c  $gfa.pangenie_genotyping.vcf.gz -t $sdf -o $gfa.pangenie_genotyping.vcf.eval
rtg vcfeval -b $truthvcf.gz -c  $gfa.pangenie_phasing.vcf.gz -t $sdf -o $gfa.pangenie_phasing.vcf.eval

echo "extract the diplotype into a pair of FASTAs (vcf2fasta)"

vcf2fasta -f $fasta.ref.fa <(zcat $gfa.pangenie_genotyping.vcf.gz)

echo "map simulated reads from pirs to both haplotypes, choosing the best match for each"
# plan: use straight mapping and ignore the mq
# todo: use a script to pick the best haplotype for each read (pair?)
echo "call variants for each haplotype"
echo "correct high-quality homozygous calls to the variant allele and emit the FASTAs"

