import os
import pandas as pd
import numpy as np
import argparse

''' argument parsing '''
print('parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
parser.add_argument('-hetdata',
                    help='path to heterogeneous data file. Only heterogeneous loci are included in the analysis.',
                    default='/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/het.txt')
parser.add_argument('-biasdata',
                    help='path to bias data file. Biased loci are excluded from the analysis.',
                    default='data/bias.txt')
parser.add_argument('-refdata',
                    help='path to reference read counts data file.',
                    default='/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ref.txt')
parser.add_argument('-altdata',
                    help='path to alternate read counts data file.',
                    default='/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/alt.txt')
parser.add_argument('-asedest',
                    help='output path to calculated ase.',
                    default='results/ase.txt')
parser.add_argument('-validdest',
                    help='output path to ase validity.',
                    default='results/ase_validity.txt')
parser.add_argument('-minreads',
                    help='minimum number of reads per ase site.',
                    type = int,
                    default=20)
parser.add_argument('-pseudo',
                    help='pseudocount',
                    type = float,
                    default=50.0)

args = parser.parse_args()


''' setting variables '''
home_dir = args.home
het_data_path = args.hetdata
bias_data_path = args.biasdata
ref_data_path = args.refdata
alt_data_path = args.altdata

ase_data_dest_path = args.asedest
validity_data_dest_path =args.validdest

MIN_READS = args.minreads
PSEUDO_COUNT = args.pseudo + 0.0

# set working directory
os.chdir(home_dir)

print('reading dgn data ...')
het_data = pd.read_table(het_data_path, sep=' ', header=None, index_col=0)
ref_data = pd.read_table(ref_data_path, sep=' ', header=None, index_col=0)
alt_data = pd.read_table(alt_data_path, sep=' ', header=None, index_col=0)
bias_data = pd.read_table(bias_data_path, sep='\t', header=None, index_col=None)

print('rearrange data to have same sample-order')
ref_data = ref_data.loc[het_data.index, :]
alt_data = alt_data.loc[het_data.index, :]

print('calculating ase ...')
delta = 0.01
total_reads = ref_data + alt_data
total_reads = total_reads.replace(0, delta) # to avoid division by zero
ase_data = abs((ref_data+PSEUDO_COUNT) / (total_reads+2.0*PSEUDO_COUNT) - 0.5)

print('saving ase data ...')
ase_data.to_csv(ase_data_dest_path, sep='\t', index=True, header=False)

print('creating ase validity matrix ...')
unbias_list = [1-row[1][0] for row in bias_data.iterrows()]
unbias_mat = np.matrix(unbias_list * ase_data.shape[0])
unbias_mat = unbias_mat.reshape(ase_data.shape[0],ase_data.shape[1])
ase_validity = het_data * (total_reads >= MIN_READS) * unbias_mat
ase_validity.to_csv(validity_data_dest_path, sep='\t', index=True, header=False)


