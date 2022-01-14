#!/usr/bin/env bash

vcftools \
--gzvcf variants.32samples.0.5missing.recombined.vcf.gz \
--chr $1 \
--recode \
--out variants.32samples.0.5missing.recombined.${1}.vcf
