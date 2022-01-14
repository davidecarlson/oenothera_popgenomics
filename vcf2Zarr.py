#!/usr/bin/env python

import allel
import numpy as np

# convert vcfs to ZARR with scikit allel

vcf = '../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz'

chrom = '../chrom_list.txt'

data = np.loadtxt(chrom, delimiter = '\n\n', dtype = str)

for chrom in data:
    allel.vcf_to_zarr(vcf, '../vcfs/variants.30samples.0.5missing.biallelicSNPs.zarr', group = chrom, region = chrom, fields = '*')
