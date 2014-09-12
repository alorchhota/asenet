import os
import pandas as pd
from scipy import stats
import numpy as np
import time

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
het_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/het.txt'
ase_data_path = 'data/ase.txt'
tf_expr_data_path = 'data/tf_expr_data.txt'


MIN_HET_SAMPLES = 30 #minimum number of samples heterogeneous at a locus

curTime = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
tfs_loci_dest_path = 'results/significant_tfs_' + curTime + '.txt'
test_statistics_dest_path = 'results/test_statistics_' + curTime + '.txt'

# set working directory
os.chdir(home_dir)

print('reading dgn ase data ...')
ase_data = pd.read_table(ase_data_path, sep='\t', header=None, index_col=0)
het_data = pd.read_table(het_data_path, sep=' ', header=None, index_col=0)

print('reading TF expression data ...')
tf_expr_data = pd.read_table(tf_expr_data_path, sep='\t', header=0, index_col=0)
tfs = tf_expr_data.columns.values


print('processing data for correlations ...')
n_loc = ase_data.shape[1]
n_tfs = tf_expr_data.shape[1]
n_sample = tf_expr_data.shape[0]

n_het_samples = [sum(het_at_locus) for _, het_at_locus in het_data.iteritems()]
n_loci_tested = sum([1 for x in n_het_samples if x >= MIN_HET_SAMPLES])

p_thresh = 0.05 / (n_loci_tested * n_tfs)

print('cleaning previous results ...')
with open(tfs_loci_dest_path, 'w') as outFile:
    outFile.write('\t'.join(['tf','locus_index','n_het_sample','r','p']) + '\n')

print('calculating correlations ...')

minLocusIndex = 1000
maxLocusIndex = 5000  # n_loc

for li in range(minLocusIndex, maxLocusIndex):
    if li % 100 == 0:
        print(li)
    
    if n_het_samples[li] < MIN_HET_SAMPLES:
        continue
    
    het_index = np.where(het_data.iloc[:,li]==1)[0]
    ase_expr = ase_data.iloc[het_index,li]
    
    # store test statistics
    test_statistics = []
    
    for tfi in range(n_tfs):
        tfg = tfs[tfi]
        tf_expr = tf_expr_data.iloc[het_index, tfi]
        
        cor = stats.spearmanr(tf_expr, ase_expr)
        test_statistics.append((tfg, li, cor[0], cor[1]))
        if cor[1] <= p_thresh:
            with open(tfs_loci_dest_path, 'a') as outFile:
                toPrint = [tfg, li, n_het_samples[li], cor[0], cor[1]]
                outFile.write('\t'.join([str(item) for item in toPrint]) + '\n')
                print(':::: ' + str(li) + ',' + tfg)

    with open(test_statistics_dest_path, 'a') as statOutFile:
        statOutFile.write( '\n'.join('\t'.join([str(item) for item in tup]) for tup in test_statistics))
        statOutFile.write( '\n')
    