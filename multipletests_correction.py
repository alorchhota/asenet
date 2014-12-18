import os
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
test_statistics_path = 'results/test_statistics_pseudo_2014-12-17_15-55-45.txt'

# method = 'fdr_bh' or 'bonferroni'
method_param = 'fdr_bh'
alpha = 0.05;

multitest_corrected_test_statistics_dest_path = 'results/test_statistics_pseudo_2014-12-17_15-55-45_' + method_param + '.txt'


# set working directory
os.chdir(home_dir)

print('reading test statistics ...')
test_stat_data = pd.read_table(test_statistics_path, sep='\t', header=0, index_col=None)

print('multiple test (' + method_param +') correction ...')
multitest_significance = multipletests(test_stat_data['p'], alpha=alpha, method=method_param)[0];
sig_indexes = np.where(multitest_significance == True)[0]

print('saving significant test statistics ...')
sig_test_stat_data = test_stat_data.iloc[sig_indexes,:]
sig_test_stat_data.to_csv(path_or_buf=multitest_corrected_test_statistics_dest_path, sep='\t', header=True, index=False)

