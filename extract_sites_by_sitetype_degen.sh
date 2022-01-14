#!/usr/bin/env bash

# split vcf into two new files based on zero-fold and four-fold codon degeneracy


vcftools \
--positions ../siteTypes.4fold.txt \
--gzvcf ../vcfs/variants.30samples.0.5missing.recombined.minDP20.vcf.gz \
--recode \
--recode-INFO-all \
--out ../vcfs/variants.30samples.0.5missing.recombined.minDP20.sitetype.4folddegen

vcftools \
--positions ../siteTypes.0fold.txt \
--gzvcf ../vcfs/variants.30samples.0.5missing.recombined.minDP20.vcf.gz \
--recode \
--recode-INFO-all \
--out ../vcfs/variants.30samples.0.5missing.recombined.minDP20.sitetype.0folddegen
