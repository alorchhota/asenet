import os
import pandas as pd
from scipy import stats
import numpy as np
import argparse

''' argument parsing '''
print('parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
parser.add_argument('-validdata',
                    help='path to validity data file. Only valid loci are included in the analysis.',
                    default='data/ase_validity.txt')
parser.add_argument('-asedata',
                    help='path to ase data file.',
                    default='data/ase.txt')
parser.add_argument('-tfdata',
                    help='path to TF expression data file.',
                    default='data/tf_expr_data.txt')
parser.add_argument('-statdest',
                    help='output path to test statistic.',
                    default='results/test_statistics.txt')
parser.add_argument('-sigdest',
                    help='output path to significant tf-asesite pairs.',
                    default='results/significant_tfs.txt')
parser.add_argument('-minsamples',
                    help='minimum number of valid samples per ase site.',
                    type = int,
                    default=30)
args = parser.parse_args()

''' setting variables '''
home_dir = args.home
validity_data_path = args.validdata
ase_data_path = args.asedata
tf_expr_data_path = args.tfdata
MIN_VALID_SAMPLES = args.minsamples
tfs_loci_dest_path = args.sigdest
test_statistics_dest_path = args.statdest

# set working directory
os.chdir(home_dir)

print('reading dgn ase data ...')
ase_data = pd.read_table(ase_data_path, sep='\t', header=None, index_col=0)
validity_data = pd.read_table(validity_data_path, sep='\t', header=None, index_col=0)

print('reading TF expression data ...')
tf_expr_data = pd.read_table(tf_expr_data_path, sep='\t', header=0, index_col=0)
tfs = tf_expr_data.columns.values


print('processing data for correlations ...')
n_loc = ase_data.shape[1]
n_tfs = tf_expr_data.shape[1]
n_sample = tf_expr_data.shape[0]

n_valid_samples = [sum(valid_at_locus) for _, valid_at_locus in validity_data.iteritems()]
n_loci_tested = sum([1 for x in n_valid_samples if x >= MIN_VALID_SAMPLES])

p_thresh = 0.05 / (n_loci_tested * n_tfs)

print('cleaning previous results ...')
with open(tfs_loci_dest_path, 'w') as outFile:
    outFile.write('\t'.join(['tf','locus_index','n_het_sample','r','p']) + '\n')
    
with open(test_statistics_dest_path, 'w') as outFile:
    outFile.write('\t'.join(['tf','locus_index','r','p']) + '\n')


print('calculating correlations ...')

minLocusIndex = 0
maxLocusIndex = n_loc

test_statistics = []

for li in range(minLocusIndex, maxLocusIndex):
    if li>0 and li % 100 == 0:
        print(li)
        with open(test_statistics_dest_path, 'a') as statOutFile:
            statOutFile.write( '\n'.join('\t'.join([str(item) for item in tup]) for tup in test_statistics))
            statOutFile.write( '\n')
        test_statistics = []
    
    if n_valid_samples[li] < MIN_VALID_SAMPLES:
        continue
    
    het_index = np.where(validity_data.iloc[:,li]==1)[0]
    ase_expr = ase_data.iloc[het_index,li]
    
    # store test statistics
    for tfi in range(n_tfs):
        tfg = tfs[tfi]
        tf_expr = tf_expr_data.iloc[het_index, tfi]
        
        cor = stats.spearmanr(tf_expr, ase_expr)
        test_statistics.append((tfg, li, cor[0], cor[1]))
        if cor[1] <= p_thresh:
            with open(tfs_loci_dest_path, 'a') as outFile:
                toPrint = [tfg, li, n_valid_samples[li], cor[0], cor[1]]
                outFile.write('\t'.join([str(item) for item in toPrint]) + '\n')
                print(':::: ' + str(li) + ',' + tfg)

if len(test_statistics) > 0:
    with open(test_statistics_dest_path, 'a') as statOutFile:
        statOutFile.write( '\n'.join('\t'.join([str(item) for item in tup]) for tup in test_statistics))
        statOutFile.write( '\n')
    

