#!/usr/bin/env python

import pandas as pd
import seaborn as sns
from scipy import stats
from matplotlib import pyplot as plt

reads = pd.read_csv('../read_counts.txt', header = None, names = ['Sample', 'Read Count'])


reads['Species'] = [i.split('_')[-1] for i in reads['Sample']]

#print(reads)

missing = pd.read_csv('../stats/variant_stats.imiss', header = 0, sep = '\t')
snps = pd.read_csv('../stats/snps_persample.txt', header = None, sep = '\t', names = ['Sample', 'SNP Count'])
#print(missing)
#print(snps)

reads['Missing'] = missing['F_MISS']
reads['SNPs'] = snps['SNP Count']
#print(reads)

means = reads.groupby(['Species']).mean()
snp_std = reads.groupby(['Species']).std()
print(means)
print(snp_std)


biennis_df = reads[reads['Species'] == 'biennis']
#print(biennis_df)
elata_df = reads[reads['Species'] == 'elata']
#print(elata_df)

snp_count_test = stats.ttest_ind(elata_df['SNPs'], biennis_df['SNPs'], permutations = 10000)

print(f'The results of 10,000 permutation tests of SNP counts between the two species are:  {snp_count_test}')

def get_corr(x, y):
    correlation = stats.pearsonr(x,y)
    return(correlation)

elata_reads_miss_cor = get_corr(elata_df['Read Count'], elata_df['Missing'])
elata_reads_snp_cor = get_corr(elata_df['Read Count'], elata_df['SNPs'])
biennis_reads_miss_cor = get_corr(biennis_df['Read Count'], biennis_df['Missing'])
biennis_reads_snp_cor = get_corr(biennis_df['Read Count'], biennis_df['SNPs'])

print(f'Correlation between read count and genotype missingness in elata is: {elata_reads_miss_cor}')
print(f'Correlation between read count and genotype missingness in biennis is: {biennis_reads_miss_cor}')

print(f'Correlation between read count and snp count in elata is: {elata_reads_snp_cor}')
print(f'Correlation between read count and snp count in biennis is: {biennis_reads_snp_cor}')
def scatterplot(df, x, y, hue, xlab, ylab,prefix):
    sns.set(rc={'figure.figsize':(4,4)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    ax = sns.scatterplot(data = df, x = x, y = y, hue = hue, palette= ['r', 'b'])
    
    #plt.xlim(0,20000000)
    #plt.ylim(0,0.6)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(prefix + '.png', bbox_inches='tight')
    plt.savefig(prefix + '.eps', bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()

scatterplot(reads, x='Read Count', y = 'Missing', hue = 'Species', xlab = 'Number of reads', ylab = 'Fraction of missing genotype calls', prefix = '../figures/reads_v_miss')

scatterplot(reads, x='Read Count', y = 'SNPs', hue = 'Species', xlab = 'Number of reads', ylab = 'Number of SNPs', prefix = '../figures/reads_v_snp_count')
