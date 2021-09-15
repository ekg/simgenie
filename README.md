# simgenie

A hold-one-out test harness for diploid haplotype inference and polishing.

Input: GFA-format graph with each haplotype stored as a path (`P`-line).

## Overview

We take an input GFA and convert it to a phased VCF against a selected reference sequence.
We remove a specific sample from this VCF, generating a haplotype panel.
Using `PanGenie`, we reconstruct the held-out sample haplotypes using this panel and reads simulated from the true diplotype.
These haplotypes will have errors caused by the inference process and the presence of rare variants in our sample.
We finally polish each estimated haplotype using freebayes, generating a consensus diplotype that we compare to our true diplotype.

## Steps

1. extract phased VCF from the GFA
2. extract selected haplotypes into FASTA
3. simulate reads (with pirs) from the diploid FASTA
4. apply pangenie to infer the underlying diplotype
5. extract the diplotype into a pair of FASTAs (vcf2fasta)
6. map simulated reads fro pirs to both haplotypes, choosing the best match for each
7. call variants for each haplotype
8. correct high-quality homozygous calls to the variant allele and emit the FASTAs

## Usage

To run the simulation using chm13 as reference and HG03579 as the target sample:

```
mkdir -p test
cd test
../gen_sim.sh ../genes/HTT.grch38.gfa.gz chm13 HG03579 HG03579#1 HG03579#2 16
```
