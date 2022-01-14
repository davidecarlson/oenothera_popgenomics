#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
import statistics
dxy =  pd.read_csv('../pixy_results/biennis_elata_pixy_500kb_dxy.txt', header = 0, sep = '\t')

#print(dxy)
#print(dxy.info())
#print(len(dxy['avg_dxy']))

pie = pd.read_csv('../pixy_results/biennis_elata_pixy_500kb_pi.txt', header = 0, sep = '\t')

#print(pie.info())

#print(pie)

biennis_pie = pie[pie['pop'] == 'biennis'].reset_index()

#biennis_pie = biennis_pie.notna()




divergence = dxy['avg_dxy']
div_mean = np.mean(divergence)
div_std = np.std(divergence)
div_med = np.nanmedian(divergence)
print(f'median divergence:      {div_med}')
print(f'standard deviation divergence:      {div_std}')
print(f'average divergence:      {div_mean}')
biennis_pie = biennis_pie.join(divergence)
#print(biennis_pie)
biennis_pie['avg_dxy'].replace('', np.nan, inplace=True)
biennis_pie.dropna(subset=['avg_dxy'], inplace=True)

biennis_pie.to_csv('test.txt', sep = '\t', index = False)

bi_corr = stats.pearsonr(biennis_pie['avg_pi'], biennis_pie['avg_dxy'])

print(f'The correlation between pi and dxy in biennis is {bi_corr[0]} with a p-value of {bi_corr[1]}')

#biennis_pie.to_csv('test.txt', sep = '\t', index = True)
#print(biennis_pie.info())


elata_pie = pie[pie['pop'] == 'elata'].reset_index()

elata_pie = elata_pie.join(divergence)
#print(elata_pie)

elata_pie['avg_dxy'].replace('', np.nan, inplace=True)
elata_pie.dropna(subset=['avg_dxy'], inplace=True)

ela_corr = stats.pearsonr(elata_pie['avg_pi'], elata_pie['avg_dxy'])
print(f'The correlation between pi and dxy in elata is {ela_corr[0]} with a p-value of {ela_corr[1]}')

combined = pd.concat([biennis_pie, elata_pie])

#print(combined)

combined.to_csv('../pixy_results/combined_pi_dxy.txt', sep ='\t', index = False)

def scatterplot(df, x, y, xlab, ylab, png):
    sns.set(rc={'figure.figsize':(6,6)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    ax = sns.scatterplot(data = df, x = x, y = y, hue = 'pop',  hue_order = ['biennis', 'elata',], style = 'pop', palette = ['r', 'b'])
    # make identity line
    mn = min(df[x].min(), df[y].min())
    mx = max(df[x].max(), df[y].max())
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(points, points, color = 'k', marker = None,
        linestyle = '--', linewidth = 1.0)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend([],[], frameon=False)
    plt.savefig('../figures/' + png + '.png', bbox_inches='tight')
    plt.savefig('../figures/' + png + '.eps', bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()


scatterplot(combined, x='avg_pi', y = 'avg_dxy', xlab = 'Nucleotide diversity (\u03C0)', ylab = 'Absolute divergence (Dxy)', png = 'pi_v_dxy')
#scatterplot(biennis_pie, x='avg_pi', y = 'avg_dxy', col = 'r', xlab = 'Nucleotide diversity (pi)', ylab = 'Absolute divergence (Dxy)', png = 'biennis_pi_v_dxy.png')
#scatterplot(elata_pie, x='avg_pi', y = 'avg_dxy', col = 'b', xlab = 'Nucleotide diversity (pi)', ylab = 'Absolute divergence (Dxy)', png = 'elata_pi_v_dxy.png')
