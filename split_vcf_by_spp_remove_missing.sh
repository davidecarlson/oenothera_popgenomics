#!/usr/bin/env bash

#split vcf file by specie and require at least one non-ref alleles

VCF=../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz

bcftools view \
--samples-file ../elata_samples.txt $VCF |\
bcftools view -g ^miss - | \
bcftools view --min-ac 1 | \
vcffixup - | bgzip -c > ../vcfs/elata.nomissing.biallelicSNPs.vcf.gz

bcftools view \
--samples-file ../biennis_samples.txt $VCF |\
bcftools view -g ^miss - | \
bcftools view --min-ac 1 | \
vcffixup - | bgzip -c > ../vcfs/biennis.nomissing.biallelicSNPs.vcf.gz
