#!/usr/bin/env bash

# Apply site filters to VCF

MISSING=0.5
QUAL=30
MINDP=5
MAXDP=100

# separate VCF into variant and invariant sites and filter separately

# get invariant sites with missingness filter applied
# added step to recalculated AC and AF


vcftools --gzvcf ../vcfs/variants.nogrand.30samples.vcf.gz \
--max-maf 0 \
--remove-filtered-all \
--max-missing $MISSING \
--recode-INFO-all \
--min-meanDP $MINDP \
--max-meanDP $MAXDP \
--recode --stdout | vcffixup - | bgzip -c > ../vcfs/variants.30samples.0.5missing.invariant.vcf.gz

tabix ../vcfs/variants.30samples.0.5missing.invariant.vcf.gz

# get variant sites and apply additional filters
vcftools --gzvcf ../vcfs/variants.nogrand.30samples.vcf.gz \
--mac 1 \
--minQ 30 \
--remove-indels \
--remove-filtered-all \
--max-missing $MISSING \
--min-meanDP $MINDP \
--max-meanDP $MAXDP \
--recode-INFO-all \
--recode --stdout | vcffixup - |bgzip -c > ../vcfs/variants.30samples.0.5missing.variant.vcf.gz

tabix ../vcfs/variants.30samples.0.5missing.variant.vcf.gz

# recombine the two VCFs

bcftools concat \
--allow-overlaps \
../vcfs/variants.30samples.0.5missing.invariant.vcf.gz \
../vcfs/variants.30samples.0.5missing.variant.vcf.gz \
-O z \
-o ../vcfs/variants.30samples.0.5missing.recombined.vcf.gz

tabix ../vcfs/variants.30samples.0.5missing.recombined.vcf.gz
