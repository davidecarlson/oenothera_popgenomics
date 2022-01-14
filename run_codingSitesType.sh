#!/usr/bin/env bash

# use the codingSitesTypes.py script from Simon Martin to get functional categories for SNPs in vcf

python /home/progs/genomics_general/codingSiteTypes.py \
-a /bowman/datahome_migrated/oenothera/annotation/genes_with_introns.gff3 \
-f gff3 \
--ignoreConflicts \
-r /bowman/datahome_migrated/oenothera/reference/elata_reference_assembly.fasta \
-o ../siteTypes.txt 2>&1 | tee ../logs/siteTypes.log
