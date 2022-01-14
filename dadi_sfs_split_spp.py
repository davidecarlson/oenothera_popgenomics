#!/usr/bin/env python

import dadi
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

dd_elata = dadi.Misc.make_data_dict_vcf("../vcfs/elata.nomissing.includesinvariant.vcf.gz", "../../popmap_nogrand.txt")
dd_4f_elata = dadi.Misc.make_data_dict_vcf("../vcfs/elata.nomissing.includesinvariant.4folddegen.recode.vcf.gz", "../../popmap_nogrand.txt")
dd_0f_elata = dadi.Misc.make_data_dict_vcf("../vcfs/elata.nomissing.includesinvariant.0folddegen.recode.vcf.gz", "../../popmap_nogrand.txt")

dd_biennis = dadi.Misc.make_data_dict_vcf("../vcfs/biennis.nomissing.includesinvariant.vcf.gz", "../../popmap_nogrand.txt")
dd_4f_biennis = dadi.Misc.make_data_dict_vcf("../vcfs/biennis.nomissing.includesinvariant.4folddegen.recode.vcf.gz", "../../popmap_nogrand.txt")
dd_0f_biennis = dadi.Misc.make_data_dict_vcf("../vcfs/biennis.nomissing.includesinvariant.0folddegen.recode.vcf.gz", "../../popmap_nogrand.txt")


def bootstrap(dd, chunksize, projection, repsize,popID):
    # take data dict and get bootstrapped  replicated of the sfs
    chunk_size = chunksize
    chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
    bootstrap_sfs = dadi.Misc.bootstraps_from_dd_chunks(chunks, repsize, 
    pop_ids = [popID], projections = [projection], mask_corners=False, polarized = False)
    return(bootstrap_sfs)

def unmask_sfs(sfs):
    # get sfs with corners unmasked
    unmasked_sfs = []
    for i in sfs:
        i.mask = False
        unmasked_sfs.append(i)
        print(i)
    return(unmasked_sfs)

def prop_sfs(sfs, limit):
    # get proportions for each allele count category and remove count categories with no data
    prop_sfs = []
    for i in sfs:
        i = i[0:limit]
        prop = [x / sum(i) for x in i]
        prop_sfs.append(prop)
    return(prop_sfs)

def make_csv(sfs, output):
    # take full unmasked bootstrapped SFS and write out csv file
    df = pd.DataFrame(sfs)
    df.to_csv(output, index = False, header = False)
    
def make_df(sfs, species, category):
    # make df of prop sfs
    df = pd.DataFrame(sfs)
    # use melt to convert from wide to long (tidy) df format
    tidy = df.melt(var_name='Allele Count', value_name= 'Proportion')
    newdf = tidy.assign(Species=species, Category=category)
    print(newdf)
    return(newdf)

def plot_sfs(dfs, cols, outpng):
    # combine input dfs
    combined = pd.concat(dfs,ignore_index=True)
    print(combined)
    
    sns.set(rc={'figure.figsize':(15,10)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    sns.lineplot(data=combined, x="Allele Count", y="Proportion", 
    ci = "sd", hue = "Species", style = "Category", palette = cols)
    plt.savefig(outpng, bbox_inches='tight')
    plt.show()
    plt.clf()


# get bootstrap replicated for each sfs
biennis_boots = bootstrap(dd_biennis, 5000000, 32,5, 'biennis')
biennis_4f_boots = bootstrap(dd_4f_biennis, 5000000,32,  5, 'biennis')
biennis_0f_boots = bootstrap(dd_0f_biennis, 5000000,32,  5, 'biennis')
elata_boots = bootstrap(dd_elata, 5000000, 32,5, 'elata')
elata_4f_boots = bootstrap(dd_4f_elata, 5000000,32, 5, 'elata')
elata_0f_boots = bootstrap(dd_0f_elata, 5000000,32, 5, 'elata')

#print(biennis_boots)
#print(elata_boots)



# unmask each set of bs reps
biennis_allsites_unmasked = unmask_sfs(biennis_boots)
biennis_4f_unmasked = unmask_sfs(biennis_4f_boots)
biennis_0f_unmasked = unmask_sfs(biennis_0f_boots)
elata_allsites_unmasked = unmask_sfs(elata_boots)
elata_4f_unmasked = unmask_sfs(elata_4f_boots)
elata_0f_unmasked = unmask_sfs(elata_0f_boots)

# convert count totals to proportions for each set of bs reps
biennis_prop_allsites = prop_sfs(biennis_allsites_unmasked, 17)
biennis_prop_4f = prop_sfs(biennis_4f_unmasked,17)
biennis_prop_0f = prop_sfs(biennis_0f_unmasked,17)
elata_prop_allsites = prop_sfs(elata_allsites_unmasked,17)
elata_prop_4f = prop_sfs(elata_4f_unmasked,17)
elata_prop_0f = prop_sfs(elata_0f_unmasked,17)

print(biennis_prop_allsites)
print(elata_prop_allsites)


# make dfs for each set of proportion bs reps
biennis_allsites_df = make_df(biennis_prop_allsites, "biennis", "All")
biennis_4f_df = make_df(biennis_prop_4f, "biennis", "4fold")
biennis_0f_df = make_df(biennis_prop_0f, "biennis", "0fold")
elata_allsites_df = make_df(elata_prop_allsites, "elata", "All")
elata_4f_df = make_df(elata_prop_4f, "elata", "4fold")
elata_0f_df = make_df(elata_prop_0f, "elata", "0fold")

# writ out a csv file for each set of count sfs reps
make_csv(biennis_allsites_unmasked, '../easysfs/biennis_allsites_sfs_bs.csv')
make_csv(biennis_4f_unmasked, '../easysfs/biennis_4f_sfs_bs.csv')
make_csv(biennis_0f_unmasked, '../easysfs/biennis_0f_sfs_bs.csv')
make_csv(elata_allsites_unmasked, '../easysfs/elata_allsites_sfs_bs.csv')
make_csv(elata_4f_unmasked, '../easysfs/elata_4f_sfs_bs.csv')
make_csv(elata_0f_unmasked, '../easysfs/elata_0f_sfs_bs.csv')

# make spectra plots

plot_sfs([biennis_allsites_df,biennis_4f_df,biennis_0f_df,elata_allsites_df,elata_4f_df,elata_0f_df], ['r', 'b'], "../figures/bootstraps_sfs.png")
plot_sfs([biennis_allsites_df,elata_allsites_df],['r', 'b'],"../figures/allsites_bootstraps_sfs.png")
plot_sfs([biennis_4f_df,biennis_0f_df,elata_4f_df, elata_0f_df],['r', 'b'],"../figures/codingsites_bootstraps_sfs.png")
plot_sfs([biennis_allsites_df,biennis_4f_df,biennis_0f_df],['r'], "../figures/biennis_bootstraps_sfs.png")
plot_sfs([elata_allsites_df,elata_4f_df,elata_0f_df], ['b'], "../figures/elata_bootstraps_sfs.png")

