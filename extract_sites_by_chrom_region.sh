#!/usr/bin/env bash

bcftools view \
--regions-file ../chrom_int_regions.tsv \
-Oz \
-o ../vcfs/variants.32samples.0.5missing.biallelicSNPs.chromint.vcf.gz \
../vcfs/variants.32samples.0.5missing.biallelicSNPs.vcf.gz

bcftools view \
--regions-file ../chrom_ext_regions.tsv \
-Oz \
-o ../vcfs/variants.32samples.0.5missing.biallelicSNPs.chromext.vcf.gz \
../vcfs/variants.32samples.0.5missing.biallelicSNPs.vcf.gz
