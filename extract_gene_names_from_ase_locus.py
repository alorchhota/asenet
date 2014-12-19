import os
import pandas as pd

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
gencode_annot_path = 'data/gencode.v19.annotation.processed.gtf'
ase_loc_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ase_locus_annot.txt'

significant_tfs_path = 'results/test_statistics_pseudo_2014-12-19_00-02-51_bonferroni.txt'
annotated_significant_tfs_dest_path = 'results/significant_tf_ase_bonferroni.txt'


# set working directory
os.chdir(home_dir)

# import files from current directory

import gtf

print('reading gencode annotations ...')
gen_annot = gtf.GencodeGTFReader(gencode_annot_path, processed=True)
gen_annot.data = gen_annot.data.loc[:, ['chr', 'feature', 'start', 'end', 'gene_name']]

print('reading dgn ase locus data ...')
ase_loc_ann_data = pd.read_table(ase_loc_path, sep=' ', header=0, index_col=None)

print('significant tfs-ase correlation ...')
sig_tf_data = pd.read_table(significant_tfs_path, sep='\t', header=0, index_col=None)

print('extracting locus gene names ...')
locus_genes = [gen_annot.get_gene_name_from_locus(
                    ase_loc_ann_data.iloc[tup['locus_index'],0], 
                    ase_loc_ann_data.iloc[tup['locus_index'],1]) 
               for index, tup in sig_tf_data.iterrows()]
 
print('saving significant tf with locus names ...')
sig_tf_data['gene'] = locus_genes
sig_tf_data.to_csv(path_or_buf=annotated_significant_tfs_dest_path, sep='\t', header=True, index=False)
