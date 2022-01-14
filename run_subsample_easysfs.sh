#!/usr/bin/env bash

# run easy_sfs on subsampled vcfs

export PATH=/home/progs/easySFS:$PATH

ls ../vcfs/subsamples/*0folddegen.vcf | parallel -j 20 --verbose "mkdir -p ../easy_sfs/subsample/{/.} && easySFS.py --proj 16,20 -a -i {} -p ../popmap_nogrand.txt -o ../easy_sfs/subsample/{/.} -v"

ls ../vcfs/subsamples/*4folddegen.vcf | parallel -j 20 --verbose "mkdir -p ../easy_sfs/subsample/{/.} && easySFS.py --proj 14,20 -a -i {} -p ../popmap_nogrand.txt -o ../easy_sfs/subsample/{/.} -v"
