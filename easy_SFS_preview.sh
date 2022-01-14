#!/usr/bin/env bash

# run easySFS in preview mode to choose projection values

export PATH=/home/progs/easySFS:$PATH

easySFS.py \
--preview \
-a \
-i $1 \
-p ../popmap_nogrand.txt \
-o $2 \
-v 2>&1 |tee ../logs/${2}_preview.log
