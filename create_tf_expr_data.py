import os
import pandas as pd

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
tf_list_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/tf_nrg2538/tf_processed.txt'
expr_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/data_used_for_eqtl_study/trans_data.txt'
sample_label_path = 'data/sample_labels.txt'

tf_expr_data_dest_path = 'results/tf_expr_data.txt'

# set working directory
os.chdir(home_dir)


print('reading expresson data ...')
expr_data = pd.read_table(expr_data_path, sep='\t', header=0, index_col=0)


print('taking only TF expressions and sort samples ...')

tfs = pd.read_table(tf_list_path, sep='\t', header=None, names=['tf']).iloc[:,0]
dgn_genes = expr_data.columns.values
tfs_in_dgn = [t for t in tfs if t in dgn_genes]
samples_table = pd.read_table(sample_label_path, sep='\t', header=None, names=['sample'])
samples = [s for s in samples_table.iloc[:,0]]
expr_data = expr_data.loc[samples, tfs_in_dgn]

print('saving expressions ...')
expr_data.to_csv(path_or_buf=tf_expr_data_dest_path, sep='\t', header=True, index=True)