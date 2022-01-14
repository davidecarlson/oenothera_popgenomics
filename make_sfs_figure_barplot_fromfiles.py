#!/usr/bin/env python

import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def read_files(input):
    # read in file and make df
    df = pd.read_csv(input, header = 0, sep = '\t')
    print(df)
    return(df)



def plot_sfs(dfs, cols, outpng):
    # combine input dfs
    combined = pd.concat(dfs,ignore_index=True)
    print(combined)
    
    sns.set(rc={'figure.figsize':(15,10)})
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.25)
    sns.catplot(data=combined, x="Allele Count", y="Proportion", 
    ci = 95, hue = "Species", col = "Category",kind="bar", palette = cols)
    plt.savefig(outpng + ".eps", bbox_inches='tight')
    plt.savefig(outpng + ".png", bbox_inches='tight')
    #plt.show()
    plt.clf()
    
def main():

    # specify input files
    biennis_4f = '../dadi/biennis_4f_equalproj.txt'
    biennis_0f = '../dadi/biennis_0f_equalproj.txt'
    biennis_all = '../dadi/biennis_allsites_equalproj.txt'
    elata_4f = '../dadi/elata_4f_equalproj.txt'
    elata_0f = '../dadi/elata_0f_equalproj.txt'
    elata_all = '../dadi/elata_allsites_equalproj.txt'
    
    # read in files
    biennis_4f_df = read_files(biennis_4f)
    biennis_0f_df = read_files(biennis_0f)
    biennis_all_df = read_files(biennis_all)
    elata_4f_df = read_files(elata_4f)
    elata_0f_df = read_files(elata_0f)
    elata_all_df = read_files(elata_all)
    
    # plot the spectra
    plot_sfs([biennis_4f_df,biennis_0f_df,elata_4f_df, elata_0f_df],['r', 'b',],"../figures/sfs_barplot_coding_equalproj")
    plot_sfs([biennis_all_df,elata_all_df],['r', 'b',],"../figures/sfs_barplot_allsites_equalproj")
    plot_sfs([biennis_4f_df,biennis_0f_df,elata_4f_df, elata_0f_df],['r', 'b',],"../figures/sfs_barplot_coding_equalproj")
    plot_sfs([biennis_all_df,elata_all_df],['r', 'b',],"../figures/sfs_barplot_allsites_equalproj")


if __name__ == "__main__":
    main()
