#!/usr/bin/env python

import allel

# convert vcfs to HDF5 with scikit allel

vcf = '../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz'

allel.vcf_to_hdf5(vcf, '../vcfs/variants.30samples.0.5missing.biallelicSNPs.h5', fields = '*')
