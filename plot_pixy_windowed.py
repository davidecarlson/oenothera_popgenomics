#!/usr/bin/env python

# plot pixy popgen sliding window results

import argparse
import sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Takes pixy results file and plots the results across a scaffold")

parser.add_argument('--input', required=True, help='Pixy results file', action='store')
parser.add_argument('--chrom', required=True, help='Which chromosome or scaffold to plot', action='store')
parser.add_argument('--prefix', required=True, help='Output image prefix', action='store')

args=parser.parse_args()

input = args.input
chrom = args.chrom
prefix = args.prefix

# read in results file, create DF, and filter out sites where there are no genotypes
def read_file(input, chrom):
    results = pd.read_csv(input, sep="\t")
    if ('_pi' in input or '_dxy' in input):
        results = results[(results['no_sites'] > 0 ) & (results['chromosome'] == chrom)]
    elif '_fst' in input:
        results = results[(results['no_snps'] > 0 ) & (results['chromosome'] == chrom)]
    else:
        print("No valid results file found")
    #results.pivot(index = 'chromosome', columns = 'pop', values = 'pop')
    return(results)


def choose_metric(input):
    if '_pi' in input:
        metric = 'avg_pi'
        label = 'Average \u03C0'
    elif '_dxy' in input:
        metric = 'avg_dxy'
        label = 'Average Dxy'
    elif '_fst' in input:
        metric = 'avg_wc_fst'
        label = 'Average Fst' 
    return(label, metric)
    
    

def line_plot(input, df, x, png):
    plot_label,plot_metric = choose_metric(input)
    print(plot_label,plot_metric)
    sns.set(rc={'figure.figsize':(8,2)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.25)
    
    if plot_metric == 'avg_dxy' or plot_metric == 'avg_wc_fst':
        sns.lineplot(data = df, x = x, y = plot_metric, color = 'purple')
        
    else:
        sns.lineplot(data = df, x = x, y = plot_metric, hue = 'pop', hue_order = ['biennis', 'elata'], palette = ['r', 'b'])
        plt.legend([],[], frameon=False)
    plt.xlabel('Scaffold Position')
    plt.ylabel(plot_label)
    plt.savefig(png + ".png", bbox_inches='tight')
    plt.savefig(png + ".eps", bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()
    
myresults = read_file(input, chrom)

print(myresults)

line_plot(input, myresults, 'window_pos_1', prefix)
	

