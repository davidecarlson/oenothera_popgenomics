#!/usr/bin/env bash

#split vcf file by specie and require at least one non-ref alleles

VCF=../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz

bcftools view \
--samples-file ../elata_samples.txt \
--min-ac 1 \
$VCF | vcffixup - | bgzip -c > ../vcfs/elata.biallelicSNPs.vcf.gz

bcftools view \
--samples-file ../biennis_samples.txt \
--min-ac 1 \
$VCF | vcffixup - | bgzip -c > ../vcfs/biennis.biallelicSNPs.vcf.gz
