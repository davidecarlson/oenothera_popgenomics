#!/usr/bin/env python

# calculate and plot Tajima's D in windows across the 5 largests scaffolds

import numpy as np
import scipy
import pandas
import seaborn as sns
import h5py
import zarr
import allel
from matplotlib import pyplot as plt
import pandas as pd

samples = samples = pd.read_csv('../popmap_nogrand.txt', sep='\t', header = None, names = ['Sample', 'Population'])
sample_selection = samples.Population.isin({'biennis', 'elata'}).values

populations = {
    'all': list(range(len(samples))),
    'biennis': samples[samples.Population == 'biennis'].index.tolist(),
    'elata': samples[samples.Population == 'elata'].index.tolist(),
}

#print(populations)


zarr_path = '../vcfs/variants.30samples.0.5missing.biallelicSNPs.zarr'


chrom = 'HiC_scaffold_16'

callset = zarr.open_group(zarr_path, mode='r')



variants = allel.VariantChunkedTable(callset[chrom]['variants'], 
                                     names=['POS', 'REF', 'ALT', 'DP', 'MQ', 'QD'],
                                     index='POS')


genotypes  = callset[chrom + '/calldata/GT']
#genotypes_subset = genotypes.subset(variants, sample_selection)
#print(genotypes.info)

gts = allel.GenotypeArray(genotypes)
#print(gts)

ac_pops = gts.count_alleles_subpops(populations, max_allele=20)
#print(ac_pops)







# plot Tajima's D

pos = variants['POS'][:]

for pop in 'all', 'biennis', 'elata':
    #print(pop, ac_pops[pop].count_segregating())
    print(f'Tajimas D over {chrom} for {pop} is: {allel.tajima_d(ac_pops[pop], pos)}')
    print(f'Wattersons theta over {chrom} for {pop} is: {allel.watterson_theta(pos=pos,ac=ac_pops[pop])}')

#range, windows, counts  = allel.windowed_statistic(pos, pos, statistic=lambda v: [v[0]], size=500000, step = 50000)

#x = np.asarray(windows).mean(axis=1)



y1, _, _ = allel.windowed_tajima_d(pos, ac_pops['biennis'][:], size = 500000, step = 500000)
y2, _, _ = allel.windowed_tajima_d(pos, ac_pops['elata'][:], size = 500000, step = 500000)

w1, _, _, _ = allel.windowed_watterson_theta(pos, ac_pops['biennis'][:], size = 500000, step = 500000)
w2, _, _, _ = allel.windowed_watterson_theta(pos, ac_pops['elata'][:], size = 500000, step = 500000)

print(f"The average value of Tajima's D on scaffold {chrom} in biennis is {np.mean(y1)} with a standard deviation of {np.std(y1)}")
print(f"The average value of Tajima's D on scaffold {chrom} in elata is {np.mean(y2)} with a standard deviation of {np.std(y2)}")

print(f"The average value of Watterson's theta on scaffold {chrom} in biennis is {np.mean(w1)} with a standard deviation of {np.std(w1)}")
print(f"The average value of Watterson's theta on scaffold {chrom} in elata is {np.mean(w2)} with a standard deviation of {np.std(w2)}")
#print(y1)
#print(y2)
'''
sns.set(rc={'figure.figsize':(8,2)})
sns.set_style("white")
sns.set_context("paper", font_scale=1.25)

#fig, ax = plt.subplots()
#sns.despine(ax=ax, top = False, right = False, bottom = False, left = False, offset = {'right' : 10, 'left': 10})
sns.lineplot(x, y1,  color = 'r')
sns.lineplot(x, y2,  color = 'b')
plt.ylabel("Tajima's $D$")
plt.xlabel('Scaffold position (bp)')
plt.xlim(0, pos.max())
#ax.legend(loc='upper left', bbox_to_anchor=(1, 1));
plt.legend([],[], frameon=False)
plt.savefig("../figures/" + chrom + "_windowed_tajima.png", bbox_inches='tight')
plt.savefig("../figures/" +chrom + "_windowed_tajima.eps", bbox_inches='tight')
plt.show()
plt.clf()
'''

