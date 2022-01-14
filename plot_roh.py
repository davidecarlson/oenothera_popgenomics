#!/usr/bin/env python

# plot pixy popgen sliding window results

import statistics as stats
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from statannot import add_stat_annotation





# read in results file, create DF, and filter out sites where there are no genotypes
def read_file(input):
    results = pd.read_csv(input, sep="\t", skiprows = [0,1,2], header = 0, 
        names = ['RG', 'Sample', 'Chrom', 'Start', 'End', 'Length', 'Number of SNPs', 'Qual'])
    results = results[(results['Qual'] >= 20 )]
    return(results)
    
elata_roh = read_file('../roh_results/elata_roh.txt')

elata_med_count = elata_roh.groupby('Sample').size().median()

elata_counts = elata_roh.groupby('Sample').size()

print(elata_counts)
print(elata_med_count)

elata_roh.insert(0, 'Species', 'elata')

elata_roh_per_samp = len(elata_roh['Length']) / 15

elata_roh_avg_len = stats.mean(elata_roh['Length'])

elata_roh_std_len = stats.stdev(elata_roh['Length'])

elata_roh_med_len = stats.median(elata_roh['Length'])

elata_roh_max_len = max(elata_roh['Length'])

elata_roh_min_len = min(elata_roh['Length'])

biennis_roh = read_file('../roh_results/biennis_roh.txt')

biennis_med_count = biennis_roh.groupby('Sample').size().median()
biennis_counts = biennis_roh.groupby('Sample').size()
print(biennis_counts)

biennis_roh.insert(0, 'Species', 'biennis')

biennis_roh_per_samp = len(biennis_roh['Length']) / 15

biennis_roh_avg_len = stats.mean(biennis_roh['Length'])

biennis_roh_std_len = stats.stdev(biennis_roh['Length'])

biennis_roh_med_len = stats.median(biennis_roh['Length'])

biennis_roh_max_len = max(biennis_roh['Length'])

biennis_roh_min_len = min(biennis_roh['Length'])

print(f'The median number of rohs per sample in elata is {elata_med_count}')
print(f'The median number of rohs per sample in biennis is {biennis_med_count}')
print(f'elata samples have {elata_roh_per_samp} runs of homozygosity on average')
print(f'biennis samples have {biennis_roh_per_samp} runs of homozygosity on average')

print(f'The average length of RoH in elata is {elata_roh_avg_len} bp with a standard deviation of {elata_roh_std_len}. The median length is {elata_roh_med_len} bp. ')
print(f'The average length of RoH in biennis is {biennis_roh_avg_len} bp with a standard deviation of {biennis_roh_std_len}. The median length is {biennis_roh_med_len} bp.')

print(f'elata RoH min: {elata_roh_min_len}')
print(f'elata RoH max: {elata_roh_max_len}')

print(f'biennis RoH min: {biennis_roh_min_len}')
print(f'biennis RoH max: {biennis_roh_max_len}')
#print(biennis_roh.head())

combined = pd.concat([biennis_roh,elata_roh])

#print(combined)
    


    
def hist_bottom(df, x, hue, xlab, ylab,prefix):
    sns.set(rc={'figure.figsize':(7,4)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.25)
    
    #f, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace':0.2})
    sns.histplot(data = df, x = x, hue = hue, palette= ['r', 'b'], stat = "probability")
    #sns.histplot(data = df, x = x, hue = hue, palette= ['r', 'b'], stat = "probability", ax = ax_top)
    #sns.histplot(data = df, x = x, hue = hue, palette= ['r', 'b'], stat = "probability", ax = ax_bottom)
    
    #sns.despine(ax=ax_bottom)
    #sns.despine(ax=ax_top, bottom=True)
    
    #ax_top.set_ylim(0.15, 0.3)
    #ax_bottom.set_ylim(0,0.07)
    
    #ax = ax_top
    
    #d = 0.02
    
    #kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    #ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    
    #ax2 = ax_bottom
    #kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    #ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        
    plt.xlim(-1000000, 10000000)
    plt.ylim(0, 0.07)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    sns.despine(top = True)
    
    #remove one of the legend
    #ax_bottom.legend_.remove()
    
    plt.savefig("../figures/" + prefix + ".png", bbox_inches='tight')
    plt.savefig("../figures/" + prefix + ".eps", bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()

def hist_top(df, x, hue, xlab, ylab,prefix):
    sns.set(rc={'figure.figsize':(7,2)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.25)
    
    #f, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace':0.2})
    sns.histplot(data = df, x = x, hue = hue, palette= ['r', 'b'], stat = "probability")

        
    plt.xlim(-1000000, 10000000)
    plt.ylim(0.2, 0.28)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    sns.despine()
    
    #remove one of the legend
    #ax_bottom.legend_.remove()
    
    plt.savefig("../figures/" + prefix + ".png", bbox_inches='tight')
    plt.savefig("../figures/" + prefix + ".eps", bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()

def violinplot(df, x, y, xlab, ylab,prefix):
    sns.set(rc={'figure.figsize':(6,3)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    ax = sns.violinplot(data = df, x = x, y = y, palette= ['r', 'b'])
    order = ['biennis', 'elata']
    add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                    box_pairs=[("biennis", "elata")],
                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig("../figures/" + prefix + ".png", bbox_inches='tight')
    plt.savefig("../figures/" + prefix + ".eps", bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()

    
#hist_bottom(combined, 'Length', 'Species', 'RoH Length (bp)', 'Proportion', 'roh_hist_bottom')
hist_top(combined, 'Length', 'Species', 'RoH Length (bp)', 'Proportion', 'roh_hist_top')    

#violinplot(combined, 'Species', 'Length', 'Species', 'RoH Length (bp)', 'roh_violin')


