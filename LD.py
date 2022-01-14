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

biennis_samples = samples.Population.isin({'biennis'}).values
elata_samples = samples.Population.isin({'elata'}).values

print(biennis_samples)
print(elata_samples)

populations = {
    'all': list(range(len(samples))),
    'biennis': samples[samples.Population == 'biennis'].index.tolist(),
    'elata': samples[samples.Population == 'elata'].index.tolist(),
}




zarr_path = '../vcfs/variants.30samples.0.5missing.biallelicSNPs.zarr'


chrom = 'HiC_scaffold_13'

callset = zarr.open_group(zarr_path, mode='r')



variants = allel.VariantChunkedTable(callset[chrom]['variants'], 
                                     names=['POS', 'REF', 'ALT', 'DP', 'MQ', 'QD'],
                                     index='POS')

# create permissive filter expression

filter_expression = '(DP > 5)'
variant_selection = variants.eval(filter_expression)[:]

genotypes  = callset[chrom + '/calldata/GT']
#genotypes_subset = genotypes.subset(variants, sample_selection)
print(genotypes.info)

gts = allel.GenotypeArray(genotypes)
print(gts)

gn = gts.to_n_alt()

print(gn)

gts_biennis = gts.subset(variant_selection, biennis_samples)
gts_elata = gts.subset(variant_selection, elata_samples)

biennis_gn = gts_biennis.to_n_alt()
elata_gn = gts_elata.to_n_alt()
    
# plot LD

pos = variants['POS'][:]

range, windows, counts  = allel.windowed_statistic(pos, pos, statistic=lambda v: [v[0]], size=1000000, step = 100000)

x = np.asarray(windows).mean(axis=1)


m = allel.rogers_huff_r(biennis_gn[:10000]) ** 2
ax = allel.plot_pairwise_ld(m)
plt.show()
plt.savefig("ld.eps", bbox_inches='tight')


y1, _, _ = allel.windowed_r_squared(pos, gn, size = 1000000, step = 100000)
print(y1)

#y2, _, _ = allel.windowed_r_squared(pos, elata_gn, size = 5000)


sns.set(rc={'figure.figsize':(8,2)})
sns.set_style("white")
sns.set_context("paper", font_scale=1.25)

#fig, ax = plt.subplots()
#sns.despine(ax=ax, top = False, right = False, bottom = False, left = False, offset = {'right' : 10, 'left': 10})
sns.lineplot(x, y1,  color = 'r')
#sns.lineplot(x, y2,  color = 'b')
plt.ylabel("Linkage disequilibrium")
plt.xlabel('Scaffold position (bp)')
plt.xlim(0, pos.max())
#ax.legend(loc='upper left', bbox_to_anchor=(1, 1));
plt.legend([],[], frameon=False)
plt.savefig("../figures/" + chrom + "_windowed_LD.png", bbox_inches='tight')
plt.savefig("../figures/" +chrom + "_windowed_LD.eps", bbox_inches='tight')
plt.show()
plt.clf()


