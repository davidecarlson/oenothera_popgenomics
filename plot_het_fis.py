#!/usr/bin/env python

# plot obs het and Fis from vcflib's popStats

import argparse
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from statannot import add_stat_annotation

parser = argparse.ArgumentParser(description="Takes popStats results file and create violin plots")

parser.add_argument('--biennis', required=True, help='biennis results file', action='store')
parser.add_argument('--elata', required=True, help='elata results file', action='store')

args=parser.parse_args()

biennis_input = args.biennis
elata_input = args.elata

# read in results file, create DF, and filter out sites where there are no genotypes
def read_files(biennis, elata):
    results_b = pd.read_csv(biennis, sep="\t", header=None, names = ['Chrom', 'Position', 'Allele Freq', 'Exp Het',
                 'Obs Het','Het Count','Hom Ref Count', 'Hom Alt Count', 'Fis'])
    
    results_b = results_b.assign(Species='biennis')
    results_e = pd.read_csv(elata, sep="\t",header=None, names = ['Chrom', 'Position', 'Allele Freq', 'Exp Het',
                 'Obs Het','Het Count','Hom Ref Count', 'Hom Alt Count', 'Fis'])
    results_e = results_e.assign(Species='elata')
    
    combined = pd.concat([results_b, results_e])

    return(combined)

myresults = read_files(biennis_input,elata_input)

myresults['Inbr Coef'] = 1 - (myresults['Obs Het'] / myresults['Exp Het'])
print(myresults)
#print(F)



def med_het(input):
    median_het = input.groupby(['Species'])[['Obs Het']].median()
    print("The median Observered Het values are:")
    print(median_het)

def mean_het(input):
    mean_het = input.groupby(['Species'])[['Obs Het']].mean()
    print("The mean Observered Het values are:")
    print(mean_het)

def std_het(input):
    std_het = input.groupby(['Species'])[['Obs Het']].std()
    print("The standard devations for Observed  Het values are:")
    print(std_het)
    
def med_f(input):
    median_f = input.groupby(['Species'])[['Inbr Coef']].median()
    print("The median Inbreeding coefficients are:")
    print(median_f)

def mean_f(input):
    mean_f = input.groupby(['Species'])[['Inbr Coef']].mean()
    print("The mean Inbreeding coefficients are:")
    print(mean_f)

def std_f(input):
    std_f = input.groupby(['Species'])[['Inbr Coef']].std()
    print("The standard devations for Inbreeding Coefficients values are:")
    print(std_f)    

med_het(myresults)
med_f(myresults)
mean_het(myresults)
mean_f(myresults)           
std_het(myresults)
std_f(myresults)

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
    plt.savefig(png + ".png", bbox_inches='tight')
    plt.savefig(png + ".eps", bbox_inches='tight')
    plt.show()
    #plt.legend(loc='right', bbox_to_anchor=(1, 0.5))
    plt.clf()
    


violinplot(myresults, 'Species', 'Obs Het', 'Species', 'Observed Heterozygosity',  '../figures/DP20_Obs_het_violin.eps')
violinplot(myresults, 'Species', 'Inbr Coef', 'Species', 'Inbreeding Coefficient (F)', '../figures/DP20_Inbreeding_coef_violing.eps')

