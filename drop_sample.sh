#!/usr/bin/env bash

# script to remove indels, filtered variants, and  drop the MTJ0790_biennis sample with extremely high missing rates
# also removing two samples that may be mislabeled

vcftools \
--gzvcf ../../../all.recal.filtered.min1geno.vcf.gz \
--remove-indels \
--remove-filtered-all \
--remove-indv MTJ0790_biennis \
--remove-indv MTJ0487_biennis \
--remove-indv OEN99_elata \
--recode \
--recode-INFO-all \
--stdout | vcffixup - |bgzip -c > ../vcfs/variants.30samples.vcf.gz
