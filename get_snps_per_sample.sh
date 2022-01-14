#!/usr/bin/env bash

FILE=../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz
for sample in $(bcftools query -l $FILE);
    do 
    count=$(bcftools view -c1 -H -s $sample $FILE | wc -l);
    printf "%s\t%s\n" "$sample" "$count";
done > ../stats/snps_persample.txt
