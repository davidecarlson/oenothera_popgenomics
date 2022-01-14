#!/usr/bin/env python

import dadi
import gzip
import os
from collections import OrderedDict
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# script to take 4fold and 0fold degenerate site VCF files, write out bootstrapped SFS files and make plots with confidence intervals

# use method of converting vcf to data dictionary modified from https://github.com/isaacovercast/easySFS to get SFS with invariant sites

def message(message):
    print(message)

def read_input(vcf_name, all_snps=True, verbose=True):

    ## Counter to track which locus we're evaluating and a list
    ## to hold all lines for each locus so we can randomly
    ## select one snp per locus if necessary
    cur_loc_number = -1
    cur_loc_snps = []

    ## use gzip? 
    read_mode = 'r'
    if vcf_name.endswith(".gz"):
        read_mode = 'rt'
        ofunc = gzip.open
    else:  
        ofunc = open
    infile = ofunc(vcf_name, read_mode)
    lines = infile.readlines()
    infile.close()

    for line in lines:
        if line.startswith("#CHROM"):
            header = line

    ## Just get the data lines, not the comments
    lines = [x for x in lines if not x.startswith('#')]
    if verbose:
        print("  Number of snps in input file: {}".format(len(lines)))

    ## lines now here has a list of either all snps in the input
    ## or a subset that includes only one snp per locus
    genotypes = pd.DataFrame([x.split() for x in lines], columns=header.split())
    return genotypes
    
def get_populations(pops_file, pop_order=[], verbose=False):
    # Here we need to read in the individual population
    # assignments file and do this:
    # - populate the locs dictionary with each incoming population name
    # - populate another dictionary with individuals assigned to populations
    # :param list pop_order: If specified, reorders the populations in the
    #   returned pops dict. This matters for dadi 3 pop analysis.
    # Add the 'U' to handle opening files in universal mode, squashes the
    # windows/mac/linux newline issue.

    try:
        with open(pops_file, 'r') as popsfile:
            ind2pop = {}
            pops = OrderedDict()
        
            lines = popsfile.readlines()
            ## Get all the populations
            for line in lines:
                pops.setdefault(line.split()[1], [])
        
            for line in lines:
                ind = line.split()[0]
                pop = line.split()[1]
                ind2pop[ind] = pop
                pops[pop].append(ind)

        print("Processing {} populations - {}".format(len( pops ), pops.keys()))
        if(verbose):
            for pop,ls in pops.items():
                print(pop, ls)

    except Exception as inst:
        msg = """
    Problem reading populations file. The file should be plain text with one
    individual name and one population name per line, separated by any amount of
    white space. There should be no header line in this file. 
    An example looks like this:
        ind1    pop1
        ind2    pop1
        ind3    pop2
        ind4    pop2"""
        print(msg)
        print("    File you specified is: ".format(pops_file))
        print("    Error - {}".format(inst))
        raise

    if pop_order:
        pop_order = [x for x in pop_order.split(",")]
        if sorted(pop_order) != sorted(list(pops.keys())):
            msg = """
    Population names in `pop_order` must be identical to those in the pops
    file. You have:
        pop_order: {}
        pops from file: {}""".format(pop_order, list(pops.keys()))
            raise Exception(msg)

        # Sort pops in order specified (py3 dicts are ordered by default).
        pops = {pop:pops[pop] for pop in pop_order}

    return ind2pop, pops

def make_datadict(genotypes, pops, verbose=False, ploidy=2):
    dd = {}

    ## Get genotype counts for each population
    for row in genotypes.iterrows():
        ## iterrows() returns a tuple for some reason
        row = row[1]

        calls = {}
        for pop in pops.keys():
            ## If there is a bunch of info associated w/ each snp then
            ## just carve it off for now.
            pop_genotypes = [row[x].split(":")[0] for x in pops[pop]]
            ref_count = sum([x == "0" or x == "0/0" or x == "0|0" for x in pop_genotypes]) * ploidy
            alt_count = sum([x == "1" or x == "1/1" or x == "1|1" for x in pop_genotypes]) * ploidy
            ## Haploids shouldn't have hets in the vcf 
            het_count = sum([x == "1/0" or x == "0/1" or x == "1|0" or x == "0|1" for x in pop_genotypes])

            ref_count += het_count
            alt_count += het_count
            calls[pop] = (ref_count, alt_count)

        dd[row["#CHROM"]+"_"+row["POS"]] =\
            {"segregating":[row["REF"], row["ALT"]],\
            "calls":calls,\
            "outgroup_allele":row["REF"]}
    return dd

# the remaining functions are ones that I wrote

def bootstrap(dd, chunksize, projection, repsize,popIDs):
    # take data dict and get bootstrapped  replicated of the sfs
    chunk_size = chunksize
    chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
    bootstrap_sfs = dadi.Misc.bootstraps_from_dd_chunks(chunks, repsize,
    pop_ids = popIDs, projections = projection, mask_corners=False, polarized = False)
    return(bootstrap_sfs)


def unmask_sfs(sfs):
    # get sfs with corners unmasked
    unmasked_sfs = []
    for i in sfs:
        i.mask = False
        unmasked_sfs.append(i)
    print(unmasked_sfs)
    return(unmasked_sfs)


def prop_sfs(sfs):
    # get proportions for each allele count category and remove count categories with no data
    prop_sfs = []
    for i in sfs:
        # get half the size of the projection plus 1
        # for this script we're keeping the invariant count category
        limit = int((len(i) + 1) / 2)
        i = i[0:limit]
        prop = [count / sum(i) for count in i]
        prop_sfs.append(prop)
    print(prop_sfs)
    return(prop_sfs)
    
    


def make_df(sfs, species, category, output):
    # make df of prop sfs
    df = pd.DataFrame(sfs)
    #df.columns = [column + 1 for column in df.columns]
    print(df)
    # use melt to convert from wide to long (tidy) df format
    tidy = df.melt(var_name='Allele Count', value_name= 'Proportion')
    newdf = tidy.assign(Species=species, Category=category)
    invariant = newdf[newdf['Allele Count'] == 0]
    invariant.to_csv(output, sep = '\t', index = False)
    print(invariant)
    return(invariant)


def plot_sfs(dfs, cols, outpng):
    # combine input dfs
    combined = pd.concat(dfs,ignore_index=True)
    #print(combined)
    
    sns.set(rc={'figure.figsize':(15,10)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    sns.catplot(data=combined, x="Allele Count", y="Proportion", 
    ci = 95, hue = "Species", col = "Category",kind="bar", palette = cols)
    plt.ylim(0.90, 1.0)
    plt.savefig(outpng + ".eps", bbox_inches='tight')
    plt.savefig(outpng + ".png", bbox_inches='tight')
    #plt.show()
    plt.clf()
   
def main():

    # specify the vcf files
    zerofold_vcf = '../vcfs/variants.30samples.0.5missing.recombined.sitetype.0folddegen.recode.vcf.gz'
    fourfold_vcf = '../vcfs/variants.30samples.0.5missing.recombined.sitetype.4folddegen.recode.vcf.gz'
    allsites_vcf = '../vcfs/variants.30samples.0.5missing.recombined.vcf.gz'
    pops_file = "../popmap_nogrand.txt"

    # convert vcf to genotypes object
    fourfold_geno = read_input(fourfold_vcf, True, True)
    zerofold_geno = read_input(zerofold_vcf, True, True)
    all_geno = read_input(allsites_vcf, True, True)


    # read in populations file
    ind2pop, pops = get_populations(pops_file,'',True)
    

    # make data dictionary from  genotypes
    dd4f = make_datadict(fourfold_geno, pops=pops, ploidy=2, verbose=True)
    dd0f = make_datadict(zerofold_geno, pops=pops, ploidy=2, verbose=True)
    ddall = make_datadict(all_geno, pops=pops, ploidy=2, verbose=True)


    # make SFS from dds. NOTE: both species projected to the same size for plotting purposes even though this is not the size that maxes Seg Sites
    elata_4f_sfs = bootstrap(dd4f, 5000000, [20], 100, ['elata'])
    elata_0f_sfs = bootstrap(dd0f, 5000000, [20], 100, ['elata'])
    biennis_4f_sfs = bootstrap(dd4f, 5000000, [14], 100, ['biennis'])
    biennis_0f_sfs = bootstrap(dd0f, 5000000, [14], 100, ['biennis'])
    elata_all_sfs = bootstrap(ddall, 5000000, [12], 100, ['elata'])
    biennis_all_sfs = bootstrap(ddall, 5000000, [18], 100, ['biennis'])
    

    # unmask the SFS

    elata_4f_unmasked = unmask_sfs(elata_4f_sfs)
    elata_0f_unmasked = unmask_sfs(elata_0f_sfs)
    biennis_4f_unmasked = unmask_sfs(biennis_4f_sfs)
    biennis_0f_unmasked = unmask_sfs(biennis_0f_sfs)
    elata_all_unmasked = unmask_sfs(elata_all_sfs)
    biennis_all_unmasked = unmask_sfs(biennis_all_sfs)
    

    # convert spectra to proprotions
    elata_4f_prop = prop_sfs(elata_4f_unmasked)
    elata_0f_prop = prop_sfs(elata_0f_unmasked)
    biennis_4f_prop = prop_sfs(biennis_4f_unmasked)
    biennis_0f_prop = prop_sfs(biennis_0f_unmasked)
    elata_all_prop = prop_sfs(elata_all_unmasked)
    biennis_all_prop = prop_sfs(biennis_all_unmasked)
    
    # make dataframes
    elata_4f_df = make_df(elata_4f_prop, "elata", "4fold", "../dadi/elata_4f_invariant.txt")
    elata_0f_df = make_df(elata_0f_prop, "elata", "0fold", "../dadi/elata_0f_invariant.txt")
    biennis_4f_df = make_df(biennis_4f_prop, "biennis", "4fold", "../dadi/biennis_4f_invariant.txt")
    biennis_0f_df = make_df(biennis_0f_prop, "biennis", "0fold", "../dadi/biennis_0f_invariant.txt")
    elata_all_df = make_df(elata_all_prop, "elata", "Allsites", "../dadi/elata_allsites_invariant.txt")
    biennis_all_df = make_df(biennis_all_prop, "biennis", "Allsites", "../dadi/biennis_allsites_invariant.txt")

    # plot the spectra
    plot_sfs([biennis_4f_df,biennis_0f_df,elata_4f_df, elata_0f_df],['r', 'b',],"../figures/sfs_barplot_coding_invariant")
    plot_sfs([biennis_all_df,elata_all_df],['r', 'b',],"../figures/sfs_barplot_allsites_invariant")

    
if __name__ == "__main__":
    main()
