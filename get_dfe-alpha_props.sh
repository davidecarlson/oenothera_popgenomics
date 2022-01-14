#!/usr/bin/env bash

# get the proportion of mutations that fall into each fitness effect category


# create directory variables to use later
PROGDIR=/home/progs/dfe-alpha-release-2.16/
DFE_ALPHA_DIR=/bowman/datahome_migrated/oenothera/genomic/gbs_2020/ddRAD/gatk/results/vcfs/drop_MTJ0790/drop2samples/dfe-alpha
DFE_RESULTS=$DFE_ALPHA_DIR/results
PROPDIR=$DFE_RESULTS/proportions
LOGDIR=$PROPDIR/logs
mkdir $PROPDIR $LOGDIR

# get proportion results in parallel

parallel \
"$PROGDIR/prop_muts_in_s_ranges \
-c $DFE_RESULTS/{1}/selected{2}_{3}epoch/est_dfe.out \
-o $PROPDIR/{1}{2}.{3}epoch_proportions.out > $LOGDIR/{1}{2}.{3}epoch_proportions.log" \
::: biennis elata ::: {1..10} ::: {1..3}


# process the results in a tab-separated summary file

for i in $PROPDIR/*out; do
    species=$(basename $i | cut -d "." -f 1| sed 's/[0-9]\+//');
    results=$(cut -d " " -f3,6,9,12 $i | tr " " "\t");
    epoch=$(basename $i | cut -d "." -f 2| cut -d "_" -f 1);
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" $species $epoch $results;
done > $PROPDIR/results_summary.txt
    
    
