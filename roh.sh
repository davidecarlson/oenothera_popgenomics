#!/usr/bin/env bash

mkdir ../roh_results

bcftools roh \
--samples-file ../elata_samples.txt \
--estimate-AF ../elata_samples.txt \
-G30 \
--skip-indels \
--output-type r \
../vcfs/variants.30samples.0.5missing.biallelicSNPs.minDP20.vcf.gz > ../roh_results/elata_roh.txt

bcftools roh \
--samples-file ../biennis_samples.txt \
--estimate-AF ../biennis_samples.txt \
-G30 \
--skip-indels \
--output-type r \
../vcfs/variants.30samples.0.5missing.biallelicSNPs.minDP20.vcf.gz > ../roh_results/biennis_roh.txt
