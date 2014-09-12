import os
import pandas as pd

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
het_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/het.txt'
sample_label_out_dest = 'results/sample_labels.txt'


# set working directory
os.chdir(home_dir)

print('reading dgn het data ...')
het_data = pd.read_table(het_data_path, sep=' ', header=None, index_col=0)

print('extracting sample labels ...')
sample_labels = het_data.index.values
with open(sample_label_out_dest, 'w') as outFile:
    outFile.write('\n'.join(sample_labels))