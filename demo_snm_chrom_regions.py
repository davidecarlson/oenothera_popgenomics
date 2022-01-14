#!/usr/bin/env python

import dadi
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

#dd = dadi.Misc.make_data_dict_vcf("../vcfs/variants.32samples.0.5missing.biallelicSNPs.vcf.gz", "../../popmap_nogrand.txt")
dd4f = dadi.Misc.make_data_dict_vcf("../vcfs/variants.32samples.0.5missing.biallelicSNPs.chromext.4folddegen.recode.vcf", "../../popmap_nogrand.txt")
dd0f = dadi.Misc.make_data_dict_vcf("../vcfs/variants.32samples.0.5missing.biallelicSNPs.chromext.0folddegen.recode.vcf", "../../popmap_nogrand.txt")

#fs_biennis = dadi.Spectrum.from_data_dict(dd, ['biennis'], projections = [18], polarized = False)
fs4f_biennis = dadi.Spectrum.from_data_dict(dd4f, ['biennis'], projections = [18], polarized = False)
fs0f_biennis = dadi.Spectrum.from_data_dict(dd0f, ['biennis'], projections = [18], polarized = False)
#fs_elata = dadi.Spectrum.from_data_dict(dd, ['elata'], projections = [20], polarized = False)
fs4f_elata = dadi.Spectrum.from_data_dict(dd4f, ['elata'], projections = [20], polarized = False)
fs0f_elata = dadi.Spectrum.from_data_dict(dd0f, ['elata'], projections = [20], polarized = False)


def get_neutral_sfs(data):
    # take in an empirical SFS and scale a standard neutral model sfs to it
    ns = data.sample_sizes 
    pts_l = [100,110,120]
    func = dadi.Demographics1D.snm
    # no need to specify any parameters since we're not optimizing anything
    params = []
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    model = func_ex(params, ns, pts_l)
    #scale model to the data
    model_scaled = dadi.Inference.optimally_scaled_sfs(model, data)
    # fold the model sfs
    model_scaled = model_scaled.fold()
    #dadi.Plotting.plot_1d_comp_multinom(model, data,fig_num=1)
    #plt.show()
    #plt.clf()
    return(model_scaled)
    
def unmask_sfs(sfs):
    # get sfs with zero allele count category unmasked and remove freq categories > 0.5
    sfs.mask[0] = False
    sfs_unmasked = [count for count in sfs if count != "--"]
    return(sfs_unmasked)

def prop_sfs(sfs):
    # get proportions for each allele count category 
    prop_sfs = [x / sum(sfs) for x in sfs]
    return(prop_sfs)


def count_alleles(sfs):
    # read in one of the sfs and get a list of allele count categories
    allele_count = [ sfs.index(i) for i in sfs]
    return(allele_count)

def make_df(snm_df, fourfold_df, zerofold_df, names):
    #take multiple frequency spectra and create a dataframe    
    allele_count = count_alleles(snm_df)
    sfs_df = pd.DataFrame(list(zip(snm_df, fourfold_df, zerofold_df, allele_count)), columns = names)
 # use melt to convert from wide to long (tidy) df format
    tidy = sfs_df.melt(id_vars=['Allele Count'],
    value_vars=['Standard Neutral', '4-fold Degenerate','0-fold Degenerate'],
    var_name = 'Spectra Category', value_name = 'Proportion')
    return(tidy)

def plot_sfs(df, species, outpng):
    if species=='biennis':
        col = 'r'
    else:
        col = 'b'
        
    sns.set(rc={'figure.figsize':(8,8)})
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.1)
    
    sns.lineplot(data=df, x="Allele Count", y="Proportion", markers=True,
     style = "Spectra Category", color = col)
    plt.title(f'Site frequency spectra plots for {species}')
    #plt.savefig(outpng, bbox_inches='tight')
    plt.show()
    plt.clf()

# get standard neutral spectra
elata_snm = get_neutral_sfs(fs4f_elata)
biennis_snm = get_neutral_sfs(fs4f_biennis)
# unmask the zero allele count category
elata_snm_unmasked = unmask_sfs(elata_snm)
elata_4f_unmasked = unmask_sfs(fs4f_elata)
elata_0f_unmasked = unmask_sfs(fs0f_elata)

biennis_snm_unmasked = unmask_sfs(biennis_snm)
biennis_4f_unmasked = unmask_sfs(fs4f_biennis)
biennis_0f_unmasked = unmask_sfs(fs0f_biennis)
# convert count sfs to proportion sfs for plotting
elata_snm_prop = prop_sfs(elata_snm_unmasked)
elata_4f_prop = prop_sfs(elata_4f_unmasked)
elata_0f_prop = prop_sfs(elata_0f_unmasked)

biennis_snm_prop = prop_sfs(biennis_snm_unmasked)
biennis_4f_prop = prop_sfs(biennis_4f_unmasked)
biennis_0f_prop = prop_sfs(biennis_0f_unmasked)

print(biennis_snm_unmasked)
print(biennis_4f_unmasked)

# make dataframe

elata_df = make_df(elata_snm_prop, elata_4f_prop,elata_0f_prop, names=['Standard Neutral', '4-fold Degenerate','0-fold Degenerate', 'Allele Count'])
biennis_df = make_df(biennis_snm_prop, biennis_4f_prop,biennis_0f_prop, names=['Standard Neutral', '4-fold Degenerate','0-fold Degenerate', 'Allele Count'])

# plot the combined standard neutral, 0-fold, and 4-fold spectra
plot_sfs(elata_df, 'elata', '../figures/elata_coding_sfs_with_neutral.png')
plot_sfs(biennis_df, 'biennis', '../figures/biennis_coding_sfs_with_neutral.png')



'''

print(model_before)

ll = dadi.Inference.ll_multinom(model, fs_elata)

theta = dadi.Inference.optimal_sfs_scaling(model, fs_elata)


p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)

#print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs_elata, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0),maxiter=50000)
#print('Finshed optimization **************************************************')

#print(f'Best-fit parameters: {popt}')

model_opt = func_ex(popt, ns, pts_l)
#print(model_opt)
model_scaled = dadi.Inference.optimally_scaled_sfs(model_opt, fs_elata)
model_scaled = model_scaled.fold()
print(model_scaled)
ll_opt = dadi.Inference.ll_multinom(model_opt, fs_elata)


theta = dadi.Inference.optimal_sfs_scaling(model_opt, fs_elata)

dadi.Plotting.plot_1d_comp_multinom(model, fs_elata,fig_num=1)
plt.show()

'''