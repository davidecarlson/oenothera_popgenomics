#!/usr/bin/env bash

# run easySFS in preview mode to choose projection values

export PATH=/home/progs/easySFS:$PATH

easySFS.py \
--proj 30,30 \
-a \
-i $1 \
-p ../popmap_nogrand.txt \
-o $2 \
-f \
-v 2>&1 | tee ${2}.log
