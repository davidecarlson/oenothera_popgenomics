#!/usr/bin/env bash

# split vcf into two new files based on zero-fold and four-fold codon degeneracy

vcftools \
--bed ../snpeff_results/snpeff_biennis_nomissing_4folddegen.bed \
--gzvcf ../vcfs/biennis.biallelicSNPs.nomissing.ann.vcf.gz \
--recode \
--recode-INFO-all \
--out ../vcfs/biennis.biallelicSNPs.nomissing.ann.4folddegen

vcftools \
--bed ../snpeff_results/snpeff_biennis_nomissing_0folddegen.bed \
--gzvcf ../vcfs/biennis.biallelicSNPs.nomissing.ann.vcf.gz \
--recode \
--recode-INFO-all \
--out ../vcfs/biennis.biallelicSNPs.nomissing.ann.0folddegen
