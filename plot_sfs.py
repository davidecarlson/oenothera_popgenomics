#!/usr/bin/env python

# plot SFS distributions made by easySFS

import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Make SFS line plots for elata and biennis")

#parser.add_argument('--indir', required=True, help='Directory input SFS files', action='store')
#parser.add_argument('--out', required=True, help='Name of output png file', action='store')
#parser.add_argument('--variant_type', required=True, help='Type of variant', action='store')
#args=parser.parse_args()
#print(args)

#dir = args.indir

# read in sfs files and create numerical list of allele frequencies

def read_sfs(file):
    # read in the file
    sfs_file = open(file)
    # get a list of the sfs elements
    counts = sfs_file.readlines()[1].strip().split(' ')[0:10]
    print(counts)
    # convert each sfs string to float
    count_float = [float(i) for i in counts]
    # get the total number of variants
    total = sum(count_float)
    # convert counts to proportions
    proportion = [count / total for count in count_float]
    print(proportion)
    return(proportion)
    

# get SFS results for both species and for all sites, 0-fold degenerate and 4-fold degenerate sites
biennis_allsites = read_sfs('../easysfs/allsites/dadi/biennis-18.sfs')

elata_allsites = read_sfs('../easysfs/allsites/dadi/elata-18.sfs')

biennis_0fold = read_sfs('../easysfs/0fold/dadi/biennis-18.sfs')
elata_0fold = read_sfs('../easysfs/0fold/dadi/elata-18.sfs')

biennis_4fold = read_sfs('../easysfs/4fold/dadi/biennis-18.sfs')
elata_4fold = read_sfs('../easysfs/4fold/dadi/elata-18.sfs')

# get the allele count categories (1 - 9 in this case)

#allele_count = [ biennis_allsites.index(i) +1 for i in biennis_allsites]
allele_count = [ biennis_allsites.index(i) for i in biennis_allsites]

sfs_df = pd.DataFrame(list(zip(biennis_allsites,elata_allsites,biennis_0fold,elata_0fold,biennis_4fold,elata_4fold,allele_count)),
 columns = ['biennis_allsites', 'elata_allsites', 'biennis_0fold', 'elata_0fold', 'biennis_4fold', 'elata_4fold', 'allele_count'])

print(sfs_df)

sns.set(rc={'figure.figsize':(9,4)})
sns.set_style("white")
sns.set_context("paper", font_scale=1.1)

ax1 = sns.lineplot(x= 'allele_count', y = 'biennis_allsites', data = sfs_df, color = 'r')
ax4 = sns.lineplot(x= 'allele_count', y = 'elata_allsites', data = sfs_df, color = 'b')

plt.xlabel('Allele Count')
plt.ylabel('Proportion')
#plt.title(f'Site frequency spectrum of {args.variant_type} variants')
plt.legend(['biennis all sites','elata all sites'])
plt.savefig('../figures/sfs_allsites.png', bbox_inches='tight')
plt.show()
plt.clf()



ax1 = sns.lineplot(x= 'allele_count', y = 'biennis_4fold', data = sfs_df, color = '#950101', linestyle="dashed", dashes=[(2, 2), (2, 2)])
ax2 = sns.lineplot(x= 'allele_count', y = 'biennis_0fold', data = sfs_df, color = '#FF5C58', linestyle="dotted", dashes=[(3, 3), (2, 2)])
ax3 = sns.lineplot(x= 'allele_count', y = 'elata_4fold', data = sfs_df, color = '#3D56B2', linestyle="dashed", dashes=[(2, 2), (2, 2)])
ax4 = sns.lineplot(x= 'allele_count', y = 'elata_0fold', data = sfs_df, color = '#5C7AEA',linestyle="dotted", dashes=[(3, 3), (2, 2)])
plt.xlabel('Allele Count')
plt.ylabel('Proportion')
#plt.title(f'Site frequency spectrum of {args.variant_type} variants')
plt.legend(['biennis 4-fold','biennis 0-fold','elata 4-fold', 'elata 0-fold'])
plt.savefig('../figures/sfs_codingsites.png', bbox_inches='tight')
plt.show()
plt.clf()