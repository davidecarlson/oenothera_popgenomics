#!/usr/bin/env bash

# run DFE-alpha on groups of previously created SFS and config files


# create directory variables to use later
PROGDIR=/home/progs/dfe-alpha-release-2.16/
DFE_ALPHA_DIR=/bowman/datahome_migrated/oenothera/genomic/gbs_2020/ddRAD/gatk/results/vcfs/drop_MTJ0790/drop2samples/dfe-alpha
DFE_CONFIGS=$DFE_ALPHA_DIR/configs
DFE_RESULTS=$DFE_ALPHA_DIR/results
BIEN_SFS_DIR=${DFE_ALPHA_DIR}/biennis_sfs
ELATA_SFS_DIR=${DFE_ALPHA_DIR}/elata_sfs
LOGDIR=${DFE_ALPHA_DIR}/logs

mkdir $LOGDIR


# run neutral single epoch analyses
parallel --verbose -j 20 "$PROGDIR/est_dfe -c $DFE_CONFIGS/{1}_neutral_group{2}_1epoch_config.txt > $LOGDIR/{1}_neutral_group{2}_1epoch.log" ::: biennis elata ::: {1..10}

# run neutral single epoch analyses
parallel --verbose -j 20 "$PROGDIR/est_dfe -c $DFE_CONFIGS/{1}_selected_group{2}_1epoch_config.txt > $LOGDIR/{1}_selected_group{2}_1epoch.log" ::: biennis elata ::: {1..10}


# run neutral 2 epoch analyses
parallel --verbose  -j 20 "$PROGDIR/est_dfe -c $DFE_CONFIGS/{1}_neutral_group{2}_2epoch_config.txt > $LOGDIR/{1}_neutral_group{2}_2epoch.log" ::: biennis elata ::: {1..10}

# run neutral 2 epoch analyses
parallel --verbose -j 20 "$PROGDIR/est_dfe -c $DFE_CONFIGS/{1}_selected_group{2}_2epoch_config.txt > $LOGDIR/{1}_selected_group{2}_2epoch.log" ::: biennis elata ::: {1..10}


# run the 3 epoch analyses in a loop because it writes out files to the same name and doing it in parallel seems to break things

#neutral
for config in $DFE_CONFIGS/*neutral*_3epoch_config.txt; do
    $PROGDIR/est_dfe -c $config;
    rm *.out;
done

# selected
for config in $DFE_CONFIGS/*selected*_3epoch_config.txt; do
    $PROGDIR/est_dfe -c $config;
    rm *.out
done
