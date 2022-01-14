#!/usr/bin/env python

# read in VCF and use Cyvcf2 to get the distribution of genotypes for each sample

from cyvcf2 import VCF
import numpy as np
import statistics as stats
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

snps = '../vcfs/variants.30samples.0.5missing.biallelicSNPs.vcf.gz'

vcf = VCF(snps, gts012 = True)



#print(samples)


def get_gts(DP):

    vcf = VCF(snps, gts012 = True)
    samples = vcf.samples

    gts = []
    for v in vcf:
        meanDP = v.INFO.get('DP') / len(samples)
        if meanDP >= DP:
            gt = np.array(v.gt_types)
            #print(gt)
            gts.append(gt)
    #print(gts)
    gt_df = pd.DataFrame(gts, columns = samples)
    #print(gt_df)
    results = gt_df.assign(minDP=DP)
    tidy = results.melt(id_vars = 'minDP', var_name = 'Sample', value_name = 'Genotype')
    return(tidy)

gts_DP5 = get_gts(5)    
gts_DP10 = get_gts(10)
gts_DP20 = get_gts(20)
gts_DP30 = get_gts(30)
gts_DP40 = get_gts(40)
gts_DP50 = get_gts(50)
gts_DP60 = get_gts(60)

#print(gts_DP5)
#print(gts_DP10)
#print(gts_DP20)


def summarize_res(df):
    samples = df['Sample'].unique()
    
    geno_summaries = []
    
    for sample in samples:
        species = sample.split("_")[1]
        minDP = df['minDP'][0]
        hom_ref_count = len(df[(df['Sample'] == sample) & (df['Genotype'] == 0)].index)
        het_count = len(df[(df['Sample'] == sample) & (df['Genotype'] == 1)].index)
        hom_var_count = len(df[(df['Sample'] == sample) & (df['Genotype'] == 2)].index)
        
        called_gts = sum([hom_ref_count, het_count, hom_var_count])
        
        hom_ref_prop = hom_ref_count / called_gts
        het_prop = het_count / called_gts
        hom_var_prop = hom_var_count / called_gts
        geno_summaries.append([sample, species, minDP, called_gts, hom_ref_count, het_count, hom_var_count, hom_ref_prop, het_prop, hom_var_prop])
    summary_df = pd.DataFrame(geno_summaries, columns = ['Sample', 'Species', 'minDP', 'GT Count',
                                'Hom Ref Count', 'Het Count', 'Hom Var Count', 'Hom Ref Prop', 'Het Prop', 'Hom Var Prop'])
    return(summary_df)
        
gt_DP5_summary = summarize_res(gts_DP5)
gt_DP10_summary = summarize_res(gts_DP10)
gt_DP20_summary = summarize_res(gts_DP20)
gt_DP30_summary = summarize_res(gts_DP30)
gt_DP40_summary = summarize_res(gts_DP40)
gt_DP50_summary = summarize_res(gts_DP50)
gt_DP60_summary = summarize_res(gts_DP60)

print(gt_DP5_summary)

def plot_gts(dfs, y, ylab, prefix):
    combined = pd.concat(dfs, ignore_index=True)
    print(combined)
    
    sns.set(rc={'figure.figsize':(12,6)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    
    ax = sns.stripplot(x="minDP", y=y, data=combined, hue = 'Species', palette = ['r', 'b'])
    plt.xlabel('Minimum Average Read Depth')
    plt.ylabel(ylab)
    plt.savefig('../figures/' + prefix + '.eps', bbox_inches='tight')
    plt.savefig('../figures/' + prefix + '.png', bbox_inches='tight')
    plt.show()
    plt.clf()

plot_gts([gt_DP5_summary,gt_DP10_summary, gt_DP20_summary, gt_DP30_summary, gt_DP40_summary, gt_DP50_summary, gt_DP60_summary]
            , "Het Prop", "Proportion of Heterozygous Genotype Calls", 'het_vs_dp')
plot_gts([gt_DP5_summary,gt_DP10_summary, gt_DP20_summary, gt_DP30_summary, gt_DP40_summary, gt_DP50_summary, gt_DP60_summary]
            , "Hom Ref Prop", "Proportion of Homozygous Reference Genotype Calls", 'homref_vs_dp')
plot_gts([gt_DP5_summary,gt_DP10_summary, gt_DP20_summary, gt_DP30_summary, gt_DP40_summary, gt_DP50_summary, gt_DP60_summary]
            , "Hom Var Prop", "Proportion of Homozygous Variant Genotype Calls", 'homvar_vs_dp')

