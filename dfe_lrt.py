#!/usr/bin/env python

# conduct LRT of nested DFE models


import glob
import pandas as pd
from scipy.stats.distributions import chi2
import statsmodels.stats.multitest as smm

biennis_results_dir = '../dfe-alpha/results/biennis/'
elata_results_dir = '../dfe-alpha/results/elata/'

biennis_1epoch = glob.glob(biennis_results_dir + 'selected*_1epoch/est_dfe.out')
biennis_2epoch = glob.glob(biennis_results_dir + 'selected*_2epoch/est_dfe.out')
biennis_3epoch = glob.glob(biennis_results_dir + 'selected*_3epoch/est_dfe.out')

elata_1epoch = glob.glob(elata_results_dir + 'selected*_1epoch/est_dfe.out')
elata_2epoch = glob.glob(elata_results_dir + 'selected*_2epoch/est_dfe.out')
elata_3epoch = glob.glob(elata_results_dir + 'selected*_3epoch/est_dfe.out')

#print(elata_3epoch)


def process_inputs(list, epochs, species):
    results = []
    for result in list:
        spp = species
        run = int(result.split('/')[4].split('_')[0].replace('selected', ''))
        with open(result) as f:
            data = f.readline().rstrip()
            lnl = float(data.split(' ')[-1])
            if epochs == 1:

                DoF = 5
                
            elif epochs == 2:

                DoF = 7
                
            elif epochs == 3:

                DoF = 9
   
            lnl_res = [spp, run, epochs, DoF, lnl]
            #print(lnl_res)
            results.append(lnl_res)
    
    #print(results)
    res_df = pd.DataFrame(results, columns = ['Species', 'Run', 'Epochs', 'DoF', 'Lnl']).sort_values('Run', ignore_index = True)
    return(res_df)


b_1ep_df = process_inputs(biennis_1epoch, 1, 'biennis')
#print(b_1ep_df)
b_2ep_df = process_inputs(biennis_2epoch, 2, 'biennis')
#print(b_2ep_df)
b_3ep_df = process_inputs(biennis_3epoch, 3, 'biennis')
#print(b_3ep_df)


e_1ep_df = process_inputs(elata_1epoch, 1, 'elata')
e_2ep_df = process_inputs(elata_2epoch, 2, 'elata')
e_3ep_df = process_inputs(elata_3epoch, 3, 'elata')

def do_LRT(df1, df2, name, species):
    LR = 2*(df2['Lnl'] - df1['Lnl'])
    dofdiff = df2['DoF'] - df1['DoF']
    d = {'Likelihood Ratio': LR, 'Degrees of Freedom': dofdiff, 'Name': name, 'Species': species}
    df = pd.DataFrame(data = d, columns = ['Likelihood Ratio', 'Degrees of Freedom', 'Name', 'Species'])
    df['pvalue'] = df.apply(lambda row: chi2.sf(row['Likelihood Ratio'], row['Degrees of Freedom']),  axis=1)
    df['padjust'] = smm.multipletests(df['pvalue'], method='b')[1]
    return(df)
 
def write_df(df, outfile):
    df.to_csv(outfile, index = False, sep = '\t')

    
b_2v1epoch = do_LRT(b_1ep_df, b_2ep_df, '2v1epoch', 'biennis')
print(b_2v1epoch)
b_3v2epoch = do_LRT(b_2ep_df, b_3ep_df, '3v2epoch', 'biennis')
print(b_3v2epoch)
b_3v1epoch = do_LRT(b_1ep_df, b_3ep_df, '3v1epoch', 'biennis')
print(b_3v1epoch)


e_2v1epoch = do_LRT(e_1ep_df, e_2ep_df, '2v1epoch', 'elata')
print(e_2v1epoch)
e_3v1epoch = do_LRT(e_1ep_df, e_3ep_df, '3v1epoch', 'elata')
print(e_3v1epoch)
e_3v2epoch = do_LRT(e_2ep_df, e_3ep_df, '3v2epoch', 'elata')


write_df(b_2v1epoch, '../dfe-alpha/results/biennis_2v1epochs.txt')
write_df(b_3v2epoch, '../dfe-alpha/results/biennis_3v2epochs.txt')
write_df(b_3v1epoch, '../dfe-alpha/results/biennis_3v1epochs.txt')
write_df(e_2v1epoch, '../dfe-alpha/results/elata_2v1epochs.txt')
write_df(e_3v2epoch, '../dfe-alpha/results/elata_3v2epochs.txt')
write_df(e_3v1epoch, '../dfe-alpha/results/elata_3v1epochs.txt')



    
