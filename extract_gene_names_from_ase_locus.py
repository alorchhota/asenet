import os
import pandas as pd
import argparse

''' argument parsing '''
print('parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
parser.add_argument('-lociannotdata',
                    help='path to ase loci annotation data file.',
                    default='data/ase_locus_gene_annot.txt')
parser.add_argument('-sigdata',
                    help='path to significant tf-asesite pairs.',
                    default='results/test_statistics_fdr_bh.txt')
parser.add_argument('-annotdest',
                    help='output path to annotated tf-asesite pairs.',
                    default='results/significant_tf_ase_fdr_bh.txt')
args = parser.parse_args()

''' setting variables '''
home_dir = args.home
ase_loc_path = args.lociannotdata
significant_tfs_path = args.sigdata
annotated_significant_tfs_dest_path = args.annotdest


# set working directory
os.chdir(home_dir)

print('reading dgn ase locus data ...')
ase_loc_ann_data = pd.read_table(ase_loc_path, sep='\t', header=0, index_col=None)

print('significant tfs-ase correlation ...')
sig_tf_data = pd.read_table(significant_tfs_path, sep='\t', header=0, index_col=None)

print('extracting locus gene names ...')
locus_genes = [ase_loc_ann_data.iloc[tup['locus_index'], 0] for index, tup in sig_tf_data.iterrows()]

print('saving significant tf with locus names ...')
sig_tf_data['gene'] = locus_genes
sig_tf_data.to_csv(path_or_buf=annotated_significant_tfs_dest_path, sep='\t', header=True, index=False)

