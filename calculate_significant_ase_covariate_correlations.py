import os
import pandas as pd
from scipy import stats
import numpy as np


''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
validity_data_path = 'data/ase_validity.txt'
ase_data_path = 'data/ase.txt'
#covariate_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/covariates/Biological_and_hidden_factors.txt'
covariate_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/covariates/Technical_factors.txt'
sample_label_path = 'data/sample_labels.txt'
sig_tf_ase_data_path = 'results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_fdr_bh.txt'


#statistic_dest_path = 'results/ase_biological_covariate_corr.txt'
statistic_dest_path = 'results/ase_technical_covariate_corr.txt'

# set working directory
os.chdir(home_dir)

print('reading dgn ase data ...')
ase_data = pd.read_table(ase_data_path, sep='\t', header=None, index_col=0)
validity_data = pd.read_table(validity_data_path, sep='\t', header=None, index_col=0)

print('reading dgn covariate data ...')
cov_data = pd.read_table(covariate_data_path, sep='\t', header=0, index_col=0)

print('reading sample label data ...')
samples_data = pd.read_table(sample_label_path, sep='\t', header=None, names=['sample'])
samples = [s for s in samples_data.iloc[:,0]]

print('reading significant tf-ase pairs data ...')
tf_ase_data = pd.read_table(sig_tf_ase_data_path, sep='\t', header=0, index_col=None)
locus_indexes = [li for li in tf_ase_data.loc[:,'locus_index']]

print('rearrange data to have same sample-order')
ase_data = ase_data.loc[samples, :]
validity_data = validity_data.loc[samples, :]
cov_data = cov_data.loc[samples, :]

print('preparing test-statistic header ...')
with open(statistic_dest_path, 'w') as outFile:
    outFile.write('\t'.join(['locus_index','gene','covariate','r','p','tf_ase_r']) + '\n')

print('calculating correlations ...')
   
# iterate each significant ase
for tf_ase_row in tf_ase_data.iterrows():
    locus_idx = tf_ase_row[1]['locus_index']
    gene = tf_ase_row[1]['gene']
    tf_ase_r = tf_ase_row[1]['r']
    # 
    # generate valid ase data
    valid_idx = np.where(validity_data.iloc[:,locus_idx]==1)[0]
    ase_vals = ase_data.iloc[valid_idx,locus_idx]
    #
    # store test statistics
    test_statistics = []
    #
    # iterate each covariate
    for cov in cov_data.columns.values:
        cov_vals = cov_data.iloc[valid_idx, :]
        cov_vals = cov_vals.loc[:, cov]
        #
        # calculate correlations and save
        cor = stats.spearmanr(ase_vals, cov_vals)
        test_statistics.append((locus_idx, gene, cov, cor[0], cor[1], tf_ase_r))    
    #
    # save correlation statistics
    with open(statistic_dest_path, 'a') as statOutFile:
        statOutFile.write( '\n'.join('\t'.join([str(item) for item in tup]) for tup in test_statistics))
        statOutFile.write( '\n')
