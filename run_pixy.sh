#!/usr/bin/env bash

# run pixy to get popgen summary stats

mkdir ../pixy_results

pixy \
--vcf ../vcfs/variants.30samples.0.5missing.recombined.vcf.gz \
--populations ../popmap_nogrand.txt \
--stats pi fst dxy \
--window_size 50000 \
--n_cores 10 \
--output_folder ../pixy_results \
--output_prefix biennis_elata_pixy_50kb 2>&1 | tee ../logs/pixy_50kb.log
