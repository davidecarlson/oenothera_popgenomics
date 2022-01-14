#!/usr/bin/env bash

# run pixy to get popgen summary stats

# get stats for each position at 1 bp resolution
#mkdir ../pixy_results

pixy \
--vcf ../vcfs/../vcfs/variants.30samples.0.5missing.recombined.vcf.gz \
--populations ../popmap_nogrand.txt \
--stats pi \
--window_size 1 \
--chunk_size 5000000 \
--n_cores 35 \
--output_folder ../pixy_results \
--output_prefix biennis_elata_pixy_allsites 2>&1 | tee ../logs/pixy_allsites.log
