''' 
This script calculates ase-covariate correlations. 
Here, every possible ASE-COV pairs are taken care of.
'''

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


#statistic_dest_path = 'results/ase_biological_covariate_corr.txt'
statistic_dest_path = 'results/ase_technical_covariate_corr.txt'

MIN_VALID_SAMPLES = 30 #minimum number of valid samples at a locus

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

print('rearrange data to have same sample-order')
ase_data = ase_data.loc[samples, :]
validity_data = validity_data.loc[samples, :]
cov_data = cov_data.loc[samples, :]

print('preparing test-statistic header ...')
with open(statistic_dest_path, 'w') as outFile:
    outFile.write('\t'.join(['locus_index','covariate','r','p']) + '\n')

print('calculating correlations ...')
# iterate each significant ase
for locus_idx in range(ase_data.shape[1]):
    if locus_idx % 1000 == 0:
        print(str(locus_idx) + ' of ' + str(ase_data.shape[1]))
    # generate valid ase data
    valid_idx = np.where(validity_data.iloc[:,locus_idx]==1)[0]
    # there must have min number of valid samples
    if len(valid_idx) < MIN_VALID_SAMPLES:
        #print('too few samples at locus idx: ' + str(locus_idx))
        continue
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
        test_statistics.append((locus_idx, cov, cor[0], cor[1]))    
    #
    # save correlation statistics
    with open(statistic_dest_path, 'a') as statOutFile:
        statOutFile.write( '\n'.join('\t'.join([str(item) for item in tup]) for tup in test_statistics))
        statOutFile.write( '\n')

