import heapq
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
tf_expr_data_path = 'data/tf_expr_data.txt'
test_statistics_path = 'results/TF-TF-Correlations-2014-09-18/tf_tf_test_statistics_2014-09-18_11-28-38_bonferroni.txt'
fig_dest_dir = 'results/TF-TF-Correlations-2014-09-18'

# set working directory
os.chdir(home_dir)

print('reading test statistics ...')
sig_test_stat_data = pd.read_table(test_statistics_path, sep='\t', header=0, index_col=None)

print('finding top n (highest r) relations')
n = 10
nthLargestR = heapq.nlargest(n,sig_test_stat_data['r'])[n-1]
sig_idx = np.where(sig_test_stat_data['r']>=nthLargestR)[0]
sig_cor = [(sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'tf1'], sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'tf2']) for idx in sig_idx]

print('reading_TF expression data ...')
tf_expr_data = pd.read_table(tf_expr_data_path, sep='\t', header=0, index_col=0)

print('generating plots ...')
for i in range(n):
    cor = sig_cor[i]
    print(cor)
    tf1_expr = tf_expr_data[cor[0]]
    tf2_expr = tf_expr_data[cor[1]]
    plt.scatter(tf1_expr, tf2_expr)
    plt.xlabel(cor[0])
    plt.ylabel(cor[1])
    plt.savefig(fig_dest_dir + "/" + str(i+1) + "_" + cor[0] + "_" + cor[1] + ".png") 
