#!/usr/bin/env bash

# interleave SFS files with mean + 95% CI intervals and reformat for use in DFE-alpha

paste -d '\n' ../easysfs/elata_0f_ci.csv ../easysfs/elata_4f_ci.csv | sed -e '1 i\3' -e '1 i\22' > ../dfe-alpha/elata_sfs.txt
paste -d '\n' ../easysfs/biennis_0f_ci.csv ../easysfs/biennis_4f_ci.csv | sed -e '1 i\3' -e '1 i\22' > ../dfe-alpha/biennis_sfs.txt
