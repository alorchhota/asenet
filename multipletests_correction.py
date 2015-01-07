import os
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
import argparse

''' argument parsing '''
print('parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
parser.add_argument('-statdata',
                    help='path to test statistic data file.',
                    default='results/test_statistics.txt')
parser.add_argument('-method',
                    help='multiple test correction method.',
                    choices=['fdr_bh','bonferroni'],
                    default='fdr_bh')
parser.add_argument('-alpha',
                    help='p-value threshold.',
                    type=float,
                    default=0.05)
parser.add_argument('-sigdest',
                    help='output path to significant tf-asesite pairs.',
                    default='results/test_statistics_fdr_bh.txt')
args = parser.parse_args()

''' setting variables '''
home_dir = args.home
test_statistics_path = args.statdata
method_param = args.method
alpha = args.alpha
multitest_corrected_test_statistics_dest_path = args.sigdest


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


