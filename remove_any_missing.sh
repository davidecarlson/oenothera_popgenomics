#!/usr/bin/env bash

# remove any sites with missing genotype calls

vcftools \
--gzvcf \
../vcfs/variants.30samples.0.5missing.recombined.vcf.gz \
--max-missing 1.0 \
--recode-INFO-all \
--recode --stdout | vcffixup - |bgzip -c > ../vcfs/variants.30samples.0.5missing.recombined.nomissing.vcf.gz
