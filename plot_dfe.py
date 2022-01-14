#!/usr/bin/env python

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats
import statsmodels.stats.multitest as smm
# read in DFE proportions summary file and calculate stats and create plots

results = '../dfe-alpha/results/proportions/results_summary.txt'

res_df = pd.read_csv(results, sep = "\t", names = ['Species', 'Epochs','Effectively Neutral', 'Mildly Deleterious', '10<NeS<100', 'NeS>100' ])

# combine the two most deleterious categories into one

combined_del = res_df['10<NeS<100'] + res_df['NeS>100']


res_df['Highly deleterious'] = combined_del

subset_df = res_df[['Species', 'Epochs','Effectively Neutral', 'Mildly Deleterious', 'Highly deleterious']]

#print(subset_df)

best = subset_df[subset_df['Epochs'] == '3epoch']
biennis_best = best[best['Species'] == 'biennis']
elata_best = best[best['Species'] == 'elata']
neutral_test = stats.mannwhitneyu(biennis_best['Effectively Neutral'], elata_best['Effectively Neutral'])
neutral_test2 = stats.ttest_ind(biennis_best['Effectively Neutral'], elata_best['Effectively Neutral'], permutations = 1000)
mild_test = stats.mannwhitneyu(biennis_best['Mildly Deleterious'], elata_best['Mildly Deleterious'])
mild_test2 = stats.ttest_ind(biennis_best['Mildly Deleterious'], elata_best['Mildly Deleterious'], permutations = 1000)
high_test = stats.mannwhitneyu(biennis_best['Highly deleterious'], elata_best['Highly deleterious'])
high_test2 = stats.ttest_ind(biennis_best['Highly deleterious'], elata_best['Highly deleterious'], permutations = 1000)

print(neutral_test2)
print(mild_test2)
print(high_test2)

pvalues = [neutral_test[1], mild_test[1], high_test[1]]

adjusted_p = smm.multipletests(pvalues, method='b')[1]

print(adjusted_p)

print(pvalues)

print(neutral_test)
print(mild_test)
print(high_test)

#print(best)

avg = subset_df.groupby(['Species', 'Epochs']).mean()


print("Average:")
print(avg)

med = subset_df.groupby(['Species', 'Epochs']).median()

print("Median:")
print(med)


print("Standard deviation:")
stdev = subset_df.groupby(['Species', 'Epochs']).std()

print(stdev)

species = subset_df.groupby(['Species']).median()

#print(species)

melted = subset_df.melt(id_vars=['Species', 'Epochs'], var_name = 'DFE Category', value_name = 'Proportion')

#print(melted)


def plot_dfe(df, cols, outpng):

    sns.set(rc={'figure.figsize':(9,6)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.25)
    sns.catplot(data=df, x="DFE Category", y="Proportion",
    hue = "Species", palette = cols, col="Epochs", sharex=True, sharey=True, s = 10)
    plt.savefig(outpng + ".png", bbox_inches='tight')
    plt.savefig(outpng + ".eps", bbox_inches='tight')
    plt.show()
    plt.clf()
    
def plot_dfe_by_spp(df, cols, outpng):

    sns.set(rc={'figure.figsize':(9,6)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    sns.catplot(data=df, x="DFE Category", y="Proportion",
    hue = "Species", palette = cols)
    plt.savefig(outpng + ".png", bbox_inches='tight')
    plt.savefig(outpng + ".eps", bbox_inches='tight')
    plt.show()
    plt.clf()
    
plot_dfe(melted, ['r', 'b'], "../figures/dfe_plot")
#plot_dfe_by_spp(melted, ['r', 'b'], "../figures/dfe_plot_by_spp.png")
