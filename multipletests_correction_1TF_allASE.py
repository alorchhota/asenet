import os
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np

''' setting variables chaning something...'''
home_dir = '/home/asaha6/github/asenet'
#home_dir = '/home/ashis/work/github/asenet'
test_statistics_path = 'results/TF-ASE-Correlation-2014-09-16/test_statistics_2014-09-15_13-06-21.txt'

# method = 'fdr_bh' or 'bonferroni'
method_param = 'bonferroni'
alpha = 0.05;

multitest_corrected_test_statistics_dest_path = 'results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_1TF_allASE_' + method_param + '.txt'


# set working directory
os.chdir(home_dir)

print('reading test statistics ...')
test_stat_data = pd.read_table(test_statistics_path, sep='\t', header=0, index_col=None)

print('grouping data by TFs...')
groups = test_stat_data.groupby(['tf']).groups
fh = open(multitest_corrected_test_statistics_dest_path, "w+")
fh.write('tf\tlocus_index\tr\tp\n');
for g in groups:
    tf_test_data = test_stat_data.iloc[groups[g],:]
    multitest_significance = multipletests(tf_test_data['p'], alpha=alpha, method=method_param)[0];
    sig_indexes = np.where(multitest_significance == True)[0]
    #
    if len(sig_indexes) > 0:
        sig_test_stat_data = tf_test_data.iloc[sig_indexes,:]
        for sig_row in sig_test_stat_data.iterrows():
            fh.write('\t'.join([str(item) for item in sig_row[1]]) + "\n")

fh.close()
