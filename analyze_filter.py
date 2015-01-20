''' this script summarizes the loci filtered by a new filter that were not filtered by the old filter '''
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
parser.add_argument('-newfdata',
                    help='path to filtered sample data by new filter.',
                    default='results/ase_validity_new.txt')
parser.add_argument('-oldfdata',
                    help='path to filtered sample data by old filter.',
                    default='results/ase_validity_old.txt')
parser.add_argument('-fildest',
                    help='output path to filtered loci.',
                    default='results/filtered_loci.txt')
parser.add_argument('-lociannotdata',
                    help='path to ase loci annotation data file.',
                    default='data/ase_locus_gene_annot.txt')

args = parser.parse_args()

''' setting variables '''
home_dir = args.home
new_data_path = args.newfdata
old_data_path = args.oldfdata
ase_loc_path = args.lociannotdata
mono_fil_dest_path = args.fildest

# set working directory
os.chdir(home_dir)

print('reading filter data ...')
new_filter_data = pd.read_table(new_data_path, sep='\t', header=None, index_col=0)
old_filter_data = pd.read_table(old_data_path, sep='\t', header=None, index_col=0)

print('reading loci annotation data ...')
ase_loc_ann_data = pd.read_table(ase_loc_path, sep='\t', header=0, index_col=None)

print('finding samples valid by old filter, but invalid by new filter ...')
fil = old_filter_data * (1-new_filter_data)

print('summarizing per loci ...')
summary = []
for li in range(new_filter_data.shape[1]):
    n_fil = sum(fil.iloc[:,li])
    if n_fil > 0:
        gene = ase_loc_ann_data.iloc[li,0]
        isPseudo = ase_loc_ann_data.iloc[li,1]
        summary.append((li, gene, isPseudo, n_fil))

print('saving loci summary ...')
with open(mono_fil_dest_path, 'w') as fh:
    fh.write('\t'.join(['loci_index', 'gene', 'is_pseudo', '#samples_filtered']) + '\n')
    text = '\n'.join(['\t'.join([str(item) for item in s]) for s in summary])
    fh.write(text)
