#!/usr/bin/env python

# plot pixy popgen sliding window results

import argparse
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from statannot import add_stat_annotation

parser = argparse.ArgumentParser(description="Takes pixy results file and create boxplot of pi values")

parser.add_argument('--input', required=True, help='Pixy results file', action='store')
parser.add_argument('--outfig', required=True, help='Output image prefix', action='store')

args=parser.parse_args()

input = args.input
outfig = args.outfig

# read in results file, create DF, and filter out sites where there are no genotypes
def read_file(input):
    results = pd.read_csv(input, sep="\t")
    results = results[(results['no_sites'] > 0 )]
    return(results)

def med_pi(input):
    median_pi = input.groupby(['pop'])[['avg_pi']].median()
    print("The median pi values are:")    
    print(median_pi)

def mean_pi(input):
    mean_pi = input.groupby(['pop'])[['avg_pi']].mean()
    print("The mean pi values are:")
    print(mean_pi)

def stdev_pi(input):
    stdev_pi = input.groupby(['pop'])[['avg_pi']].std()
    print("The STDEV pi values are:")
    print(stdev_pi)

def violinplot(df, x, y, xlab, ylab,png):
    sns.set(rc={'figure.figsize':(7,3)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    order = ['biennis', 'elata']
    ax = sns.violinplot(data = df, x = x, y = y, order = order, palette = ['r', 'b'])
    add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                    box_pairs=[("biennis", "elata")],
                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(png + '.png', bbox_inches='tight')
    plt.savefig(png + '.eps', bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()
    
myresults = read_file(input)

print(myresults.head())

med_pi(myresults)
mean_pi(myresults)
stdev_pi(myresults)

violinplot(myresults, 'pop', 'avg_pi', 'Species', 'Nucleotide Diversity (\u03C0)', outfig)
	
