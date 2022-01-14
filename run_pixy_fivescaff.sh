#!/usr/bin/env bash

# run pixy to get popgen summary stats
# providing bed file to get overlapping windows for just the 5 largest scaffolds

mkdir pixy_results

pixy \
--vcf ../vcfs/variants.30samples.0.5missing.recombined.vcf.gz \
--populations ../popmap_nogrand.txt \
--stats pi fst dxy \
--bed_file ../../../windows.bed \
--n_cores 10 \
--output_folder ../pixy_results \
--output_prefix fivescaff_500kbwindow_50kstep 2>&1 | tee logs/pixy_5scaff.log
