#!/usr/bin/env bash

# use snpEff output to get list of 4-fold degenerate sites

#modified from https://github.com/Thatguy027/SFS_Invariant_Sites/blob/master/2020_SFS_Analysis/Paper_Files/Divergent/FileXX_make_files_Divergent.sh

bcftools query -f '%CHROM\t%POS\t%END\t%ANN\n'  $1 |\
grep "synonymous_variant" |\
grep "protein_coding" |\
grep -v "splice_" |\
cut -f-4 |\
cut -f1,11,13 -d"|" |\
sed 's/[A-Z]|//g' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $0, $5 % 3}' |\
awk '{if(($4 ~ "p.Ser" || $4 ~ "p.Pro" || $4 ~ "p.Thr" || $4 ~ "p.Ala" || $4 ~ "p.Val" || $4 ~ "p.Leu" || $4 ~ "p.Gly" || $4 ~ "p.Arg") && $7 == 0) print $1, $2=$2-1, $3, "4FOLD"}' OFS="\t"
