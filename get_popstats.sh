#!/usr/bin/env bash

#version=1.0.2

# user vcflib's popStats to get Fis estimates and observed Het

popStats \
--type PL \
--target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
--file ../vcfs/elata.biallelicSNPs.vcf.gz > ../stats/elata.biallelicSNPs.popstats

popStats \
--type PL \
--target 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14 \
--file ../vcfs/biennis.biallelicSNPs.vcf.gz > ../stats/biennis.biallelicSNPs.popstats
