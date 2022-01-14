#!/usr/bin/env python

import allel
import argparse
import seaborn as sns
import numpy as np
import pandas as pd
from matplotlib  import pyplot as plt

parser = argparse.ArgumentParser(description="Get PCA plots from VCF")

parser.add_argument('--vcf', required=True, help='input VCF file', action='store')
parser.add_argument('--out_png', required=True, help='output png file', action='store')
parser.add_argument('--popmap', required=True, help='Tab separate file listing each sample and the species/pop it belongs to', action='store')

args=parser.parse_args()


# use scikit-allel to produce PCA from multi-sample VCF file

vcf_file = args.vcf
popmap = args.popmap
png = args.out_png
eps = png.replace('png', 'eps')

# make DF from popmap info
df_samples = pd.read_csv(popmap, sep = '\t', header = None)
df_samples.columns = ['Sample', 'Species']


# drop sample MTJ0562_grandiflora because it's major outlier
drop = ['MTJ0562_grandiflora','MTJ0564_grandiflora','MTJ0551_grandiflora']

df_samples_subset = df_samples[df_samples['Sample'].isin(drop) == False]

#print(removed)
#print(df_samples_subset)

mysamples = df_samples_subset['Sample'].tolist()

print(mysamples)

species = df_samples_subset.Species.unique()


# read in the data

callset = allel.read_vcf(vcf_file, fields='*', samples = mysamples)

#print(sorted(callset.keys()))

# create genotype array

gt = allel.GenotypeArray(callset['calldata/GT'])

#print(gt.shape)



# get allele count

ac = gt.count_alleles()


# remove singletons and multiallelic variants
print(np.count_nonzero(ac.max_allele() > 1))
print(np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1)))

flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
gf = gt.compress(flt, axis=0)


#transform into a 2-D array
gn = gf.to_n_alt()



# create pca without filtering for LD
coords1, model1 = allel.pca(gn, n_components=2, scaler='patterson')

print(coords1)
print(model1)



# assign colors to species

spp_colours = {
	'elata': 'b',
	'biennis': 'r',
	'grandiflora': '#008000',
	}

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for spp in species:
        flt = (sample_population == spp)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=spp_colours[spp], 
                label=spp, markersize=6, mec='k', mew=.5)
        for i, txt in enumerate(mysamples):
            ax.annotate(txt, (x[i], y[i]))
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))


def fig_pca(coords, model, title, sample_population=None):
    sns.set(rc={'figure.figsize':(4,4)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    if sample_population is None:
        sample_population = df_samples_subset.Species.values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure()
    #fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 1, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    #ax = fig.add_subplot(1, 2, 2)
    #plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    #ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()


pca_ax = fig_pca(coords1, model1, 'Conventional SNP PCA without LD Pruning')
plt.savefig(png, bbox_inches='tight')
plt.savefig(eps, bbox_inches='tight')
'''# create PCA without variance scaling
coords2, model2 = allel.pca(gn, n_components=10, scaler=None)

pca_ax2 = fig_pca(coords2, model2, 'Conventional SNP PCA without LD Pruning and variance scaling')
plt.savefig("snp_pca_noprune_noscaling.png", bbox_inches='tight')


# function for LD pruning

def ld_prune(gn, size, step, threshold=.1, n_iter=3):
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn


# try removing linked snps
gnu = ld_prune(gn, size=500, step=200, threshold=.1, n_iter=3)


# make PCA after LD pruning
coords3, model3 = allel.pca(gnu, n_components=10, scaler=None)
pca_ax3 = fig_pca(coords3, model3, 'Conventional SNP PCA with LD pruning')
plt.savefig("snp_pca_LDprune.png", bbox_inches='tight')

#print(coords1)'''
