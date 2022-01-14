#!/usr/bin/env bash

vcftools \
--gzvcf ../vcfs/variants.30samples.vcf.gz \
--remove ../grand_individuals.txt \
--recode \
--recode-INFO-all \
--stdout | vcffixup -| bgzip -c > ../vcfs/variants.nogrand.30samples.vcf.gz
