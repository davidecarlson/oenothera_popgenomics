#!/usr/bin/env bash
 
bcftools view \
--regions HiC_scaffold_2,HiC_scaffold_7,HiC_scaffold_9,HiC_scaffold_13,HiC_scaffold_16 \
-O z \
-i 'GT[*]="alt"' \
--output variants.32samples.0.5missing.recombined.5scaffs.vcf.gz \
variants.32samples.0.5missing.recombined.vcf.gz
