#!/usr/bin/env bash

# subsample each 4fold and 0fold VCF file to generate bootstrap replicates for SFS and DFE analysis

mkdir ../vcfs/subsamples

for i in {1..200};  do
    gatk SelectVariants \
    -V ../vcfs/variants.32samples.0.5missing.recombined.sitetype.4folddegen.recode.vcf.gz \
    --select-random-fraction 0.2 \
    -O ../vcfs/subsamples/subsample_${i}.sitetype.4folddegen.vcf
done
    
for i in {1..200};  do
    gatk SelectVariants \
    -V ../vcfs/variants.32samples.0.5missing.recombined.sitetype.0folddegen.recode.vcf.gz \
    --select-random-fraction 0.2 \
    -O ../vcfs/subsamples/subsample_${i}.sitetype.0folddegen.vcf
done
