import os
import pandas as pd
from scipy import stats
import numpy as np
import time

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
tf_expr_data_path = 'data/tf_expr_data.txt'

curTime = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
test_statistics_dest_path = 'results/tf_tf_test_statistics_' + curTime + '.txt'

# set working directory
os.chdir(home_dir)

print('reading TF expression data ...')
tf_expr_data = pd.read_table(tf_expr_data_path, sep='\t', header=0, index_col=0)
tfs = tf_expr_data.columns.values

print('processing data for correlations ...')
n_tfs = tf_expr_data.shape[1]

print('cleaning previous results ...')
with open(test_statistics_dest_path, 'w') as outFile:
    outFile.write('\t'.join(['tf1','tf2','r','p']) + '\n')


print('calculating correlations ...')

for tf1i in range(n_tfs-1):
    tf1 = tfs[tf1i]
    tf1_expr = tf_expr_data.iloc[:,tf1i]

    #print(tf1_expr[:5])
    #print(tf1_expr[-5:])
    
    # store test statistics
    test_statistics = []
    
    for tf2i in range(tf1i+1, n_tfs):
        tf2 = tfs[tf2i]
        tf2_expr = tf_expr_data.iloc[:, tf2i]

        #print(tf2_expr[:5])
        #print(tf2_expr[-5:])
 
        cor = stats.spearmanr(tf1_expr, tf2_expr)
        test_statistics.append((tf1, tf2, cor[0], cor[1]))

    with open(test_statistics_dest_path, 'a') as statOutFile:
        statOutFile.write( '\n'.join('\t'.join([str(item) for item in tup]) for tup in test_statistics))
        statOutFile.write( '\n')
    
