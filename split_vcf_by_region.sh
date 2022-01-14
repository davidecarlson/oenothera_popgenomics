#!/usr/bin/env bash

# select variants in 3 Mb chunks from the largest 16 scaffolds and save them into separate vcf files

cat ../16scaffolds.8mbWindows.bed | sed -e 's/\t/:/' -e 's/\t/-/' |\
while read line; do
bcftools view -O z -o ../vcfs/subsamples/$line.sitetype.0folddegen.vcf.gz \
../vcfs/variants.32samples.0.5missing.recombined.sitetype.0folddegen.recode.vcf.gz $line;
bcftools view -O z -o ../vcfs/subsamples/$line.sitetype.4folddegen.vcf.gz \
../vcfs/variants.32samples.0.5missing.recombined.sitetype.4folddegen.recode.vcf.gz $line;
done
