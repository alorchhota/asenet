import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
validity_data_path = 'data/ase_validity.txt'
ase_data_path = 'data/ase.txt'
pseudo_data_path = 'data/ase_locus_gene_annot.txt'
sample_label_path = 'data/sample_labels.txt'
MIN_VALID_SAMPLES = 30

#statistic_dest_path = 'results/ase_biological_covariate_corr.txt'
dest_fig_path = 'results/ase_distribution_pseudo.png'

# set working directory
os.chdir(home_dir)

print('reading dgn ase data ...')
ase_data = pd.read_table(ase_data_path, sep='\t', header=None, index_col=0)
validity_data = pd.read_table(validity_data_path, sep='\t', header=None, index_col=0)

print('reading sample label data ...')
samples_data = pd.read_table(sample_label_path, sep='\t', header=None, names=['sample'])
samples = [s for s in samples_data.iloc[:,0]]

print('reading pseudogene data ...')
pseudo_data = pd.read_table(pseudo_data_path, sep='\t', header=0, index_col=None)

print('rearranging data to have same sample-order ...')
ase_data = ase_data.loc[samples, :]
validity_data = validity_data.loc[samples, :]


print('taking valid ase ...')
valid_regular_ase_vals = []
valid_pseudo_ase_vals = []
valid_regular_loci_indexes = []
valid_pseudo_loci_indexes = []
for locus_idx in range(ase_data.shape[1]):
    # generate valid ase data
    if locus_idx % 1000 == 0:
        print(str(locus_idx) + ' of ' + str(ase_data.shape[1]))
    valid_idx = np.where(validity_data.iloc[:,locus_idx]==1)[0]
    # each locus must have at least min number of valid samples
    if len(valid_idx) < MIN_VALID_SAMPLES:
        continue
    #
    ase_vals = ase_data.iloc[valid_idx,locus_idx]
    if pseudo_data.iloc[locus_idx, 1] == 0:
        valid_regular_loci_indexes.append(locus_idx)
        for item in list(ase_vals):
            valid_regular_ase_vals.append(item)
    else:
        valid_pseudo_loci_indexes.append(locus_idx)
        for item in list(ase_vals):
            valid_pseudo_ase_vals.append(item)
        #valid_pseudo_ase_vals = valid_pseudo_ase_vals + list(ase_vals)

print('#valid regular loci: ' + str(len(valid_regular_loci_indexes)))
print('#valid pseudo loci: ' + str(len(valid_pseudo_loci_indexes)))

print('creating histogram plots ...')
f, subplots = plt.subplots(2,1)
# plot regular ASEs
gi=0
subplots[gi].hist(valid_regular_ase_vals, bins=100)
subplots[gi].set_title('Regular ASE Distribution')
# plot pseudo ASE
gi=1
subplots[gi].hist(valid_pseudo_ase_vals, bins=100)
subplots[gi].set_title('Pseudo ASE Distribution')
subplots[gi].set_xlabel('ase')
# save
plt.savefig(dest_fig_path)
plt.close()
