#!/usr/bin/env bash

# run pixy to get popgen summary stats

# get stats for each position at 1 bp resolution
#mkdir ../pixy_results

pixy \
--vcf ../vcfs/variants.30samples.0.5missing.recombined.sitetype.4folddegen.recode.vcf.gz \
--populations ../popmap_nogrand.txt \
--stats pi \
--window_size 1 \
--sites_file ../siteTypes.4fold.txt \
--chunk_size 5000000 \
--n_cores 35 \
--output_folder ../pixy_results \
--output_prefix biennis_elata_pixy_sitetype_4fold 2>&1 | tee ../logs/pixy_sitetype4fold.log
