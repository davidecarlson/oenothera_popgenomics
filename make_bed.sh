#!/usr/bin/env bash

mkdir ../plink_files

plink \
--vcf ../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz \
--make-bed \
--allow-extra-chr \
--snps-only \
--out ../plink_files/biallelicSNPs.bed 2>&1 | tee make_bed.log
