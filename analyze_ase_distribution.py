import os
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
validity_data_path = 'data/ase_validity.txt'
ase_data_path = 'data/ase.txt'
sample_label_path = 'data/sample_labels.txt'

#statistic_dest_path = 'results/ase_biological_covariate_corr.txt'
dest_fig_path = 'results/ase_distribution.png'

# set working directory
os.chdir(home_dir)

print('reading dgn ase data ...')
ase_data = pd.read_table(ase_data_path, sep='\t', header=None, index_col=0)
validity_data = pd.read_table(validity_data_path, sep='\t', header=None, index_col=0)

print('reading sample label data ...')
samples_data = pd.read_table(sample_label_path, sep='\t', header=None, names=['sample'])
samples = [s for s in samples_data.iloc[:,0]]

print('rearranging data to have same sample-order ...')
ase_data = ase_data.loc[samples, :]
validity_data = validity_data.loc[samples, :]

print('taking valid ase ...')
valid_ase_vals = [];
for locus_idx in range(ase_data.shape[1]):
    # generate valid ase data
    if locus_idx % 1000 == 0:
        print(str(locus_idx) + ' of ' + ase_data.shape[1])
    valid_idx = np.where(validity_data.iloc[:,locus_idx]==1)[0]
    ase_vals = ase_data.iloc[valid_idx,locus_idx]
    valid_ase_vals = valid_ase_vals + list(ase_vals);

print('creating histogram plots ...');
h = plt.hist(valid_ase_vals, bins=100)
plt.title('ASE Distribution')
plt.xlabel('ase')
plt.savefig(dest_fig_path)
plt.close()
