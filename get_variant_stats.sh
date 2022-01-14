#!/usr/bin/env bash

# get stats about the variants to inform filtering

VCF=../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz
OUT=../stats/variant_stats


#get freq of each variant
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2

# get mean depth per individual
vcftools --gzvcf $VCF --depth --out $OUT

# get mean depth per site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT

# get site quality scores for each site 
vcftools --gzvcf $VCF --site-quality --out $OUT

# get the proportion of missing data for each individual
vcftools --gzvcf $VCF --missing-indv --out $OUT

# get the proportion of missing data for each site
vcftools --gzvcf $VCF --missing-site --out $OUT

# get het and inbreeding stats
vcftools --gzvcf $VCF --het --out $OUT
