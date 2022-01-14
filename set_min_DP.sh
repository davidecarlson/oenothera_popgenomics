#!/usr/bin/env bash

minDP=20

vcftools \
--gzvcf ../vcfs/variants.30samples.0.5missing.recombined.vcf.gz \
--min-meanDP $minDP \
--recode-INFO-all \
--recode --stdout | vcffixup - |bgzip -c > ../vcfs/variants.30samples.0.5missing.recombined.minDP${minDP}.vcf.gz

tabix ../vcfs/variants.30samples.0.5missing.recombined.minDP${minDP}.vcf.gz
bcftools index ../vcfs/variants.30samples.0.5missing.recombined.minDP${minDP}.vcf.gz
