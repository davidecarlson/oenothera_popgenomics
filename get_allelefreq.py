#!/usr/bin/env python

# read in VCF and use Cyvcf2 to get the distribution of variant allele frequencies

from cyvcf2 import VCF
import numpy as np
import statistics as stats
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

#snps = '../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz'


#vcf = VCF(snps, gts012 = True)


#print(samples)


def get_freq(variants):
    vcf = VCF(variants, gts012 = True)
   
    samples = vcf.samples

    freqs = []
    for v in vcf:
 
            alt_freq = np.array(v.aaf)
            #print(alt_freq)
            freqs.append(alt_freq)
    mean_freq = np.mean(freqs)
    median_freq = np.median(freqs)
    std_freq = np.std(freqs)

    return(freqs, mean_freq, median_freq, std_freq)

biennis_freqs, biennis_mean, biennis_median, biennis_std = get_freq('../vcfs/biennis.nomissing.biallelicSNPs.vcf.gz')
print(f'The mean and median allele frequency in biennis is {biennis_mean} and {biennis_median} with a standard deviation of {biennis_std}')

elata_freqs, elata_mean, elata_median, elata_std = get_freq('../vcfs/elata.nomissing.biallelicSNPs.vcf.gz')
print(f'The mean and median allele frequency in elata is {elata_mean} and {elata_median} with a standard deviation of {elata_std}')

def get_freq_nofixed(variants):
    vcf = VCF(variants, gts012 = True)
   
    samples = vcf.samples

    freqs = []
    for v in vcf:
 
            alt_freq = np.array(v.aaf)
            #print(alt_freq)
            freqs.append(alt_freq)
    unfixed = [freq for freq in freqs if freq < 1.0]
    fixed = [freq for freq in freqs if freq == 1.0]
    #print(unfixed)
    mean_freq = np.mean(unfixed)
    median_freq = np.median(unfixed)
    std_freq = np.std(unfixed)

    return(freqs, unfixed, fixed, mean_freq, median_freq, std_freq)

biennis_freqs, biennis_unfixed, biennis_fixed, biennis_mean, biennis_median, biennis_std = get_freq_nofixed('../vcfs/biennis.nomissing.biallelicSNPs.vcf.gz')
elata_freqs, elata_unfixed, elata_fixed, elata_mean, elata_median, elata_std = get_freq_nofixed('../vcfs/elata.nomissing.biallelicSNPs.vcf.gz')
print(f'The mean and median allele frequency among VARIABLE SITES in biennis is {biennis_mean} and {biennis_median} with a standard deviation of {biennis_std}')
print(f'Biennis has {len(biennis_unfixed)} variable sites and {len(biennis_fixed)} fixed variants')
print(f'The mean and median allele frequency among VARIABLE SITES in elata is {elata_mean} and {elata_median} with a standard deviation of {elata_std}')
print(f'Elata has {len(elata_unfixed)} variable sites and {len(elata_fixed)} fixed variants')