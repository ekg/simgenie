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
bgzip -f $truthvcf && tabix -p vcf $truthvcf.gz

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
bgzip -f $gfa.pangenie_genotyping.vcf && tabix -p vcf $gfa.pangenie_genotyping.vcf.gz
bgzip -f $gfa.pangenie_phasing.vcf && tabix -p vcf $gfa.pangenie_phasing.vcf.gz


echo "use rtg eval to compare the truth and called VCFs"

sdf=$fasta.ref.sdf
rm -rf $sdf
rtg format -o $sdf $fasta.ref.fa
rm -rf $gfa.pangenie_genotyping.vcf.eval
rm -rf $gfa.pangenie_phasing.vcf.eval
rtg vcfeval -b $truthvcf.gz -c  $gfa.pangenie_genotyping.vcf.gz -t $sdf -o $gfa.pangenie_genotyping.vcf.eval
rtg vcfeval -b $truthvcf.gz -c  $gfa.pangenie_phasing.vcf.gz -t $sdf -o $gfa.pangenie_phasing.vcf.eval

echo "extract the diplotype into a pair of FASTAs (vcf2fasta)"

vcf2fasta -f $fasta.ref.fa <(zcat $gfa.pangenie_phasing.vcf.gz )

echo "map simulated reads from pirs to both haplotypes, choosing the best match for each"

hap1ref=sample_*0.fa
hap2ref=sample_*1.fa

<$hap1ref tr : _ >x && mv x $hap1ref
<$hap2ref tr : _ >x && mv x $hap2ref

samtools faidx $hap1ref
samtools faidx $hap2ref

bwa index $hap1ref
bwa index $hap2ref

bwa mem -t $threads $hap1ref $simfq | samtools sort -n >$simfq.$hap1.bynames.bam
bwa mem -t $threads $hap2ref $simfq | samtools sort -n >$simfq.$hap2.bynames.bam

# write the SAM headers
samtools view -H $simfq.$hap1.bynames.bam >$simfq.$hap1.sam
samtools view -H $simfq.$hap2.bynames.bam >$simfq.$hap2.sam

# magic best-alignment-picker script
paste -d '\n' <(samtools view $simfq.$hap1.bynames.bam) <(samtools view $simfq.$hap2.bynames.bam) \
    | awk -v seed=73 'BEGIN { srand(seed); }
                      NR > 1 && NR % 2 == 0 { if (lastq > $5) {
                         print "a", prev;
                        } else if ($5 > lastq) { print "b", $0
                        } else if (rand() > 0.5) { print "a", prev
                        } else { print "b", $0 } }
                      { prev = $0; lastq=$5; }' > $simfq.$hap1+$hap2.mix.sam

grep '^a' $simfq.$hap1+$hap2.mix.sam | tr ' ' '\t' | cut -f 2- >>$simfq.$hap1.sam
grep '^b' $simfq.$hap1+$hap2.mix.sam | tr ' ' '\t' | cut -f 2- >>$simfq.$hap2.sam

samtools view -b $simfq.$hap1.sam | samtools sort >$simfq.$hap1.bam && samtools index $simfq.$hap1.bam
samtools view -b $simfq.$hap2.sam | samtools sort >$simfq.$hap2.bam && samtools index $simfq.$hap2.bam

echo "call variants for each haplotype"

freebayes -p 2 -f $hap1ref $simfq.$hap1.bam | bcftools view --no-version -Ou >$simfq.$hap1.bcf
freebayes -p 2 -f $hap2ref $simfq.$hap2.bam | bcftools view --no-version -Ou >$simfq.$hap2.bcf

echo "correct high-quality homozygous calls to the variant allele and emit the FASTAs"

# https://github.com/VGP/vgp-assembly/blob/master/pipeline/freebayes-polish/consensus.sh

bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=$threads $simfq.$hap1.bcf >$simfq.$hap1.fb.changes.vcf.gz
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=$threads $simfq.$hap2.bcf >$simfq.$hap2.fb.changes.vcf.gz
tabix -p vcf $simfq.$hap1.fb.changes.vcf.gz
tabix -p vcf $simfq.$hap2.fb.changes.vcf.gz

bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $hap1ref $simfq.$hap1.fb.changes.vcf.gz >$simfq.$hap1.fa
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $hap2ref $simfq.$hap2.fb.changes.vcf.gz >$simfq.$hap2.fa


echo "building a reference + hap graph with pggb"

hap1name=$(head -1 $simfq.$hap1.fa | tr -d '>')
hap2name=$(head -1 $simfq.$hap2.fa | tr -d '>')

cat $fasta.ref.fa $simfq.$hap1.fa $simfq.$hap2.fa | sed "s/>$refprefix#1/>$refprefix/" | sed s/$hap1name/$sample#1#a/ | sed s/$hap2name/$sample#2#b/ >$refprefix+$sample.fa
samtools faidx $refprefix+$sample.fa

pggb -i $refprefix+$sample.fa -t $threads -p 98 -s 100000 -n 2 -k 311 -G 13117,13219 -O 0.03 -j 0 -v -U -S -o $sample.pggb

pigz -f -p $threads $sample.pggb/*.smooth.gfa

resultgfagz=$sample.pggb/*.smooth.gfa.gz

../vcf_extract_phased.sh $resultgfagz $refprefix $threads

vcfkeepsamples $refprefix*.smooth.vcf $sample | sed s/$refprefix\_// | vcffixup - | vcffilter -f 'AC > 0' | bgzip >$sample.calls.pggb.vcf.gz
tabix -p vcf $sample.calls.pggb.vcf.gz

rm -rf $sample.calls.pggb.vcf.gz.eval
rtg vcfeval -b $truthvcf.gz -c $sample.calls.pggb.vcf.gz -t $sdf -o $sample.calls.pggb.vcf.gz.eval
