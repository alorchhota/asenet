''' 
This script counts number of positive and negative tf-ase relations for every significant covariates.
'''


import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as p
from scipy import stats




''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_covariate_sig_corr_fdr_bh.txt'
frequent_cov_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_cov_fdr_sig_corr_frequent_cov.txt'
tf_ase_data_path = 'results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_fdr_bh.txt'

# output file
dest_path = 'results/cov_pos_neg_r_count.txt'

# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
cov_ase_corr_data = pd.read_table(cov_ase_corr_data_path, sep='\t', header=0, index_col=False)
frequent_cov_data = pd.read_table(frequent_cov_data_path, sep='\t', header=None, index_col=0)
tf_ase_data = pd.read_table(tf_ase_data_path, sep='\t', header=0, index_col=False)

# merge biological and technical data
#corr_data = pd.concat([bio_corr_data, tech_corr_data])

print('counting positive and negative tf-ase correlation per covariate')
tf_ase_data['tf_ase_pos'] = tf_ase_data['r'] >= 0
groups = tf_ase_data.groupby(['locus_index', 'tf_ase_pos']).groups
groupKeys = groups.keys()

counts = {}
for row in cov_ase_corr_data.iterrows():
    ase_locus = row[1]['locus_index']
    cov = row[1]['covariate']
    pos = len(groups[(ase_locus, True)]) if (ase_locus,True) in groupKeys else 0
    neg = len(groups[(ase_locus, False)]) if (ase_locus,False) in groupKeys else 0 
    if cov in counts:
        counts[cov]['pos'] += pos
        counts[cov]['neg'] += neg
    else:
       counts[cov] = {'pos': pos, 'neg':neg}


 
print(counts)

with open(dest_path, 'w') as outFile:
    for i in range(frequent_cov_data.shape[0]):
        fc = frequent_cov_data.index.values[i]
        posCount = counts[fc]['pos']
        negCount = counts[fc]['neg']
        outFile.write('\t'.join([str(fc), str(posCount), str(negCount)]) + '\n');
    



unbiased_r = []
for row in  tf_ase_data.iterrows():
    ase_locus = row[1]['locus_index']
    if ase_locus not in cov_ase_corr_data['locus_index']:
        unbiased_r.append(row[1]['r'])















