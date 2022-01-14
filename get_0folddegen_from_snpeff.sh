#!/usr/bin/env bash

# use snpEff output to get list of 4-fold degenerate sites

#modified from https://github.com/Thatguy027/SFS_Invariant_Sites/blob/master/2020_SFS_Analysis/Paper_Files/Divergent/FileXX_make_files_Divergent.sh

bcftools query -f '%CHROM\t%POS\t%END\t%ANN\n' $1 |\
grep "missense" |\
awk -F"||," '{print $1,$2,$3,$4}' OFS="\t" |\
cut -f-4 |\
cut -f1,11,13 -d"|" |\
sed 's/[A-Z]|//g' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $0, $5 % 3}' OFS="\t" |\
awk '$4 !~ "p.Ser" {print}' |\
awk '{if ($4 ~ "p.Met" || $4 ~ "p.Trp")  print $0, "0FOLD";
      else if(($4 ~ "p.Phe" || $4 ~ "p.Ile" || $4 ~ "p.Val" || $4 ~ "p.Pro" || $4 ~ "p.Thr" || $4 ~ "p.Ala" || $4 ~ "p.Tyr" || $4 ~ "p.His" || $4 ~ "p.Gln" || $4 ~ "p.Asn" || $4 ~ "p.Lys" || $4 ~ "p.Asp" || $4 ~ "p.Glu" || $4 ~ "p.Cys" || $4 ~ "p.Gly") && ($7 == 1 || $7 == 2)) print $0, "0FOLD";
      else if(($4 ~ "p.Leu" || $4 ~ "p.Arg") && $7 == 1) print $0, "0FOLD"}' OFS="\t" |\
cut -f1,2,3,8 |\
awk '{print $1, $2 = $2 - 1, $3, $4}' OFS="\t"
