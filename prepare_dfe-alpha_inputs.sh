#!/usr/bin/env bash

# use bootstrapped spectra to run DFE-alpha

# create directory variables to use later
DADI_DIR=/bowman/datahome_migrated/oenothera/genomic/gbs_2020/ddRAD/gatk/results/vcfs/drop_MTJ0790/drop2samples/dadi
DFE_ALPHA_DIR=/bowman/datahome_migrated/oenothera/genomic/gbs_2020/ddRAD/gatk/results/vcfs/drop_MTJ0790/drop2samples/dfe-alpha
DFE_CONFIGS=$DFE_ALPHA_DIR/configs
DFE_RESULTS=$DFE_ALPHA_DIR/results
SYM_SFS=$DFE_ALPHA_DIR/sfs_from_dadi
BIEN_SFS_DIR=${DFE_ALPHA_DIR}/biennis_sfs
ELATA_SFS_DIR=${DFE_ALPHA_DIR}/elata_sfs

# symlink dadi SFS files into a single directory

rm $SYM_SFS/*.sfs
find $DADI_DIR -name "*.sfs" -exec ln -s {} $SYM_SFS  ';'

# concatenate the bootstrapped spectra from each of the selected and neutral SFS files

for num in {1..100}; do

    printf "%s\n%s\n" "$(sed "2q;d" $SYM_SFS/biennis_0f_${num}.sfs)" "$(sed '2q;d' $SYM_SFS/biennis_4f_${num}.sfs)" > $SYM_SFS/biennis_combined_${num}.sfs;
    printf "%s\n%s\n" "$(sed "2q;d" $SYM_SFS/elata_0f_${num}.sfs)" "$(sed '2q;d' $SYM_SFS/elata_4f_${num}.sfs)" > $SYM_SFS/elata_combined_${num}.sfs;
done

# take combined bootstrapped spectra files and convert them into input SFS files for DFE-alpha
# each input will include 10 properly formated selected and neutral spectra

parallel -N10 "cat {} | sed -e '1 i\10' -e '1~2 i\14' > $BIEN_SFS_DIR/biennis_sfs{#}.txt " ::: $SYM_SFS/biennis_combined*
parallel -N10 "cat {} | sed -e '1 i\10' -e '1~2 i\20' > $ELATA_SFS_DIR/elata_sfs{#}.txt " ::: $SYM_SFS/elata_combined*


# make neutral config files for each run of DFE-alpha


parallel "cat > $DFE_CONFIGS/biennis_neutral_group{1}_{2}epoch_config.txt << EOF
sfs_input_file $BIEN_SFS_DIR/biennis_sfs{1}.txt
est_dfe_results_dir $DFE_RESULTS/biennis/neutral{1}_{2}epoch
data_path_2 $DFE_ALPHA_DIR/biennis_3epoch
fold 1
site_class 0
epochs {2}
search_n2 1
t2_variable 1
t2 50
EOF
" ::: {1..10} ::: {1..3}

parallel "cat > $DFE_CONFIGS/elata_neutral_group{1}_{2}epoch_config.txt << EOF
sfs_input_file $ELATA_SFS_DIR/elata_sfs{1}.txt
est_dfe_results_dir $DFE_RESULTS/elata/neutral{1}_{2}epoch
data_path_2 $DFE_ALPHA_DIR/elata_3epoch
fold 1
site_class 0
epochs {2}
search_n2 1
t2_variable 1
t2 50
EOF
" ::: {1..10} ::: {1..3}

# make selected config files for each run of DFE-alpha

parallel "cat > $DFE_CONFIGS/biennis_selected_group{1}_{2}epoch_config.txt << EOF
sfs_input_file $BIEN_SFS_DIR/biennis_sfs{1}.txt
data_path_2 $DFE_ALPHA_DIR/biennis_3epoch
est_dfe_results_dir $DFE_RESULTS/biennis/selected{1}_{2}epoch
est_dfe_demography_results_file $DFE_RESULTS/biennis/neutral{1}_{2}epoch/est_dfe.out
fold 1
epochs {2}
site_class 1
mean_s_variable 1
mean_s -0.1
beta_variable 1
beta 0.5
EOF
" ::: {1..10} ::: {1..3}

parallel "cat > $DFE_CONFIGS/elata_selected_group{1}_{2}epoch_config.txt << EOF
sfs_input_file $ELATA_SFS_DIR/elata_sfs{1}.txt
data_path_2 $DFE_ALPHA_DIR/elata_3epoch
est_dfe_results_dir $DFE_RESULTS/elata/selected{1}_{2}epoch
est_dfe_demography_results_file $DFE_RESULTS/elata/neutral{1}_{2}epoch/est_dfe.out
fold 1
epochs {2}
site_class 1
mean_s_variable 1
mean_s -0.1
beta_variable 1
beta 0.5
EOF
" ::: {1..10} ::: {1..3}

