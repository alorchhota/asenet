import os
import pandas as pd

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
het_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/het.txt'
ref_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ref.txt'
alt_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/alt.txt'

ase_data_dest_path = 'results/ase_noabs.txt'
validity_data_dest_path ='results/ase_validity_noabs.txt'

MIN_READS = 20 

# set working directory
os.chdir(home_dir)

print('reading dgn data ...')
het_data = pd.read_table(het_data_path, sep=' ', header=None, index_col=0)
ref_data = pd.read_table(ref_data_path, sep=' ', header=None, index_col=0)
alt_data = pd.read_table(alt_data_path, sep=' ', header=None, index_col=0)

print('rearrange data to have same sample-order')
ref_data = ref_data.loc[het_data.index, :]
alt_data = alt_data.loc[het_data.index, :]

print('calculating ase ...')
delta = 0.01
total_reads = ref_data + alt_data
total_reads = total_reads.replace(0, delta) # to avoid division by zero
ase_data = ref_data / total_reads - 0.5

print('saving ase data ...')
ase_data.to_csv(ase_data_dest_path, sep='\t', index=True, header=False)

print('creating ase validity matrix ...')
ase_validity = het_data * (total_reads >= MIN_READS)
ase_validity.to_csv(validity_data_dest_path, sep='\t', index=True, header=False)

