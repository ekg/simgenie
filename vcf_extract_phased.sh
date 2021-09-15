#!/bin/bash

set -eo pipefail

input=$1
ref=$2
threads=$3

prefix=$(basename $input .gfa.gz)

gfa=$(basename $input .gz)
# fixme: hack; add "haplotype identifiers" to the reference paths
# so that all paths can be loaded into the GBWT with sample annotations
zcat $input | sed 's/chm13#/chm13#1#/g' | sed 's/grch38#/grch38#1#/g' >$gfa

echo "building the GBWTGraph index for $gfa"
gbwt=$prefix.gbwt
gg=$prefix.gg
TEMPDIR=$(pwd) time vg gbwt --num-threads $threads --num-jobs $threads -G -g $gg -o $gbwt --path-regex '(.+?)#(.+?)#(.+?)' --path-fields 'CSH' $gfa

echo "building the XG index for $gfa"
xg=$prefix.xg
TEMPDIR=$(pwd) time vg convert $gg -b $gbwt -t $threads -x >$xg

echo "building VCF from $gfa and $gbwt"
vcf=$prefix.vcf
TEMPDIR=$(pwd) time vg deconstruct -a -P $ref -g $gbwt -t $threads $xg >$vcf

#rm $gfa
#pigz -p $threads $vcf
#pigz -p $threads $gbwt
#pigz -p $threads $gg
#pigz -p $threads $xg
