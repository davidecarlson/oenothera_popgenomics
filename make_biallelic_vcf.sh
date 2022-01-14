#!/usr/bin/env bash

# take filtered vcf with both invariant and variant sites and select only bialleic SNP sites

bcftools view \
-m2 \
-M2 \
-c1 \
-v snps \
-O z \
-o ../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz \
../vcfs/variants.30samples.0.5missing.recombined.vcf.gz
