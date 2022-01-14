#!/usr/bin/env bash

# split vcf into two new files based on zero-fold and four-fold codon degeneracy


vcftools \
--positions ../../elata_4f_degen.txt \
--gzvcf ../vcfs/variants.32samples.0.5missing.recombined.vcf.gz \
--recode \
--recode-INFO-all \
--out ../vcfs/variants.32samples.0.5missing.recombined.4folddegen

vcftools \
--positions ../../elata_0f_degen.txt \
--gzvcf ../vcfs/variants.32samples.0.5missing.recombined.vcf.gz \
--recode \
--recode-INFO-all \
--out ../vcfs/variants.32samples.0.5missing.recombined.0folddegen
