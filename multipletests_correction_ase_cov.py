'''
This script finds significantly correlated ase-covariate pairs.
It uses either Bonferroni or FDR correction. 
'''

import os
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
bio_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_biological_covariate_corr.txt'
tech_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_technical_covariate_corr.txt'

# method = 'fdr_bh' or 'bonferroni'
method_param = 'fdr_bh'
alpha = 0.05;

multitest_corrected_corr_dest_path = 'results/ase_covariate_sig_corr_' + method_param + '.txt'


# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
bio_corr_data = pd.read_table(bio_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)
tech_corr_data = pd.read_table(tech_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)
corr_data = bio_corr_data.append(tech_corr_data)

print('removing nan p-values')
non_nan_idx = np.where(np.isnan(corr_data['p'])==False)[0]
corr_data = corr_data.iloc[non_nan_idx,:]

print('multiple test (' + method_param +') correction ...')
multitest_significance = multipletests(corr_data['p'], alpha=alpha, method=method_param)[0];
sig_indexes = np.where(multitest_significance == True)[0]

print('saving significant test statistics ...')
sig_corr_data = corr_data.iloc[sig_indexes,:]
sig_corr_data.to_csv(path_or_buf=multitest_corrected_corr_dest_path, sep='\t', header=True, index=False)

