#!/usr/bin/env bash

# use snpEff to annotate VCF

#conda activate snp

mkdir ../snpeff_results

snpEff \
ann \
-v \
-dataDir /datahome/public/snpEff_databases \
elata_ref \
-lof \
-s ../snpeff_results/${1%.vcf.gz}.snpeff.html \
$1 > ../vcfs/${1%.vcf.gz}.ann.vcf 2> ../logs/snpeff.log
