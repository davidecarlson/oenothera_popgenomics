#!/usr/bin/env python

# get stats and comparisons of pi for synonymous and nonysynonymous sitse

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import unicodedata
from matplotlib import pyplot as plt


syn = pd.read_csv('../pixy_results/biennis_elata_pixy_sitetype_4fold_pi.noNA.txt', 
    header = 0, sep = "\t").assign(Type = 'Synonymous')
nonsyn = pd.read_csv('../pixy_results/biennis_elata_pixy_sitetype_0fold_pi.noNA.txt', 
    header = 0, sep = "\t").assign(Type = 'Nonsynonymous')

print(nonsyn[nonsyn['avg_pi'] > 0.0])

bi_syn = syn[syn["pop"] == "biennis"]
bi_nonsyn = nonsyn[nonsyn["pop"] == "biennis"]

bi_syn_zero = bi_syn[bi_syn['avg_pi'] == 0.0]
bi_syn_greater = bi_syn[bi_syn['avg_pi'] > 0.0]

bi_nonsyn_zero = bi_nonsyn[bi_nonsyn['avg_pi'] == 0.0]
bi_nonsyn_greater = bi_nonsyn[bi_nonsyn['avg_pi'] > 0.0]



elata_syn = syn[syn["pop"] == "elata"]
elata_nonsyn = nonsyn[nonsyn["pop"] == "elata"]

elata_syn_zero = elata_syn[elata_syn['avg_pi'] == 0.0]
elata_syn_greater = elata_syn[elata_syn['avg_pi'] > 0.0]

elata_nonsyn_zero = elata_nonsyn[elata_nonsyn['avg_pi'] == 0.0]
elata_nonsyn_greater = elata_nonsyn[elata_nonsyn['avg_pi'] > 0.0]

biennis_syn_avg = np.mean(bi_syn['avg_pi'])
elata_syn_avg = np.mean(elata_syn['avg_pi'])


biennis_nonsyn_avg = np.mean(bi_nonsyn['avg_pi'])
elata_nonsyn_avg = np.mean(elata_nonsyn['avg_pi'])

biennis_nonsyn_stdev = np.std(bi_nonsyn['avg_pi'])
elata_nonsyn_stdev = np.std(elata_nonsyn['avg_pi'])

biennis_syn_stdev = np.std(bi_syn['avg_pi'])
elata_syn_stdev = np.std(elata_syn['avg_pi'])

print(f'biennis mean piN:   {biennis_nonsyn_avg}')
print(f'elata mean piN:   {elata_nonsyn_avg}')
print(f'biennis mean piS:   {biennis_syn_avg}')
print(f'elata mean piS:   {elata_syn_avg}')
print(f'biennis piN standard deviation:   {biennis_nonsyn_stdev}')
print(f'elata piN standard deviation:   {elata_nonsyn_stdev}')
print(f'biennis piS standard deviation:   {biennis_syn_stdev}')
print(f'elata piS standard deviation:   {elata_syn_stdev}')

print(f'biennis piN/piS: {biennis_nonsyn_avg / biennis_syn_avg}')
print(f'elata piN/piS: {elata_nonsyn_avg / elata_syn_avg}')

print(f'In biennis, there are {len(bi_syn.index)} total syonymous sites, with {len(bi_syn_zero.index)} synonymous sites with pi = 0.0 and {len(bi_syn_greater.index)} sites with pi > 0')
print(f'In biennis, there are {len(bi_nonsyn.index)} total nonsyonymous sites, with {len(bi_nonsyn_zero.index)} nonsynonymous sites with pi = 0.0 and {len(bi_nonsyn_greater.index)} sites with pi > 0')
print(f'In elata, there are {len(elata_syn.index)} total syonymous sites, with {len(elata_syn_zero.index)} synonymous sites with pi = 0.0 and {len(elata_syn_greater.index)} sites with pi > 0')
print(f'In elata, there are {len(elata_nonsyn.index)} total nonsyonymous sites, with {len(elata_nonsyn_zero.index)} nonsynonymous sites with pi = 0.0 and {len(elata_nonsyn_greater.index)} sites with pi > 0')

syn_test = stats.mannwhitneyu(bi_syn['avg_pi'], elata_syn['avg_pi'])
print(f'comparison of synonymous sites:        {syn_test}')

nonsyn_test = stats.mannwhitneyu(bi_nonsyn['avg_pi'], elata_nonsyn['avg_pi'])
print(f'comparison of non-synonymous sites:    {nonsyn_test}')

def plot_pi(df, prefix):
    #combined = pd.concat(dfs, ignore_index=True)
    #print(combined)
    
    sns.set(rc={'figure.figsize':(12,6)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    
    ax = sns.displot(x="avg_pi", data=df, hue = 'pop', palette = ['r', 'b'], kind="ecdf")
    plt.xlabel(f'Nucleotide diversity ({unicodedata.lookup("GREEK SMALL LETTER PI")})')
    plt.xlim(-0.001, 0.001)
    plt.savefig('../figures/' + prefix + '.eps', bbox_inches='tight')
    plt.savefig('../figures/' + prefix + '.png', bbox_inches='tight')
    plt.show()
    plt.clf()

plot_pi(syn, '../figures/piN')
plot_pi(syn, '../figures/piS')
    
