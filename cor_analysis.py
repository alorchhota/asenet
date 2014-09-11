import os
import pandas as pd
from scipy import stats
import numpy as np
#import matplotlib.pyplot as plt

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
tf_data_path = 'data/tf-nrg2538-s3.txt'
gencode_annot_path = 'data/gencode.v19.annotation.processed.gtf'
het_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/het.txt'
ref_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ref.txt'
alt_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/alt.txt'
ase_loc_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ase_locus_annot.txt'
expr_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/data_used_for_eqtl_study/trans_data.txt'


MIN_HET_SAMPLES = 30 #minimum number of samples heterogeneous at a locus

out_tfs_loci = 'results/significant_tfs.txt'

# set working directory
os.chdir(home_dir)

# import files from current directory
import tfdb
import gtf

# read TFs
print('reading transcription factors ...')
tf = tfdb.TFDB(tf_data_path)

print('reading gencode annotations ...')
gen_annot = gtf.GencodeGTFReader(gencode_annot_path, processed=True)

print('reading dgn data ...')
het_data = pd.read_table(het_data_path, sep=' ', header=None, index_col=0)
ref_data = pd.read_table(ref_data_path, sep=' ', header=None, index_col=0)
alt_data = pd.read_table(alt_data_path, sep=' ', header=None, index_col=0)
expr_data = pd.read_table(expr_data_path, sep='\t', header=0, index_col=0)
ase_loc_ann_data = pd.read_table(ase_loc_path, sep=' ', header=0, index_col=None)

print('calculating ase ...')
ase_data = abs(ref_data / (ref_data + alt_data) - 0.5)
ase_data = ase_data.replace(np.Inf, 0)
ase_data = ase_data.replace(-np.Inf, 0)

print('finding genes available in dgn data ...')
gtex_genes = list(expr_data.columns.values)
tfs_in_gtex = [t for t in tf.sym if t in gtex_genes]
tf_expr_data = expr_data[tfs_in_gtex]

''' clear a few variables so that garbage collection can free memory '''
ref_data = None
alt_data = None
expr_data = None

print('calculating correlations ...')
''' remove this line'''
tfs_in_gtex = tfs_in_gtex[0:3] 

p_thresh = 0.05 / (len(tfs_in_gtex) * ase_data.shape[1])
tf_loc_pairs = []

for tfg in tfs_in_gtex:
    print(tfg)
    tf_expr = tf_expr_data[tfg]
    for li in range(ase_data.shape[1]):
        n_het_samples = sum(het_data.iloc[:,li])
        if n_het_samples < MIN_HET_SAMPLES:
            continue
        tf_ase_expr = [(e,a) for h,e,a 
                       in zip(het_data.loc[het_data.index, het_data.columns[li]], 
                              tf_expr_data.loc[het_data.index, tfg],
                              ase_data.loc[het_data.index,ase_data.columns[li]]) 
                       if h==1]
        cor = stats.spearmanr(tf_ase_expr)
        if cor[1] < p_thresh:
            locus_gene = gen_annot.get_gene_name_from_locus(ase_loc_ann_data.iloc[li,0], ase_loc_ann_data.iloc[li,1])
            tf_loc_pairs.append([tfg, locus_gene, li, n_het_samples, cor[0], cor[1]])
            print(':::: ' + locus_gene)
        


''' save significant tf and associated locations'''
print('saving significant correlations ...')
tf_df = pd.DataFrame(tf_loc_pairs, columns=['tf', 'gene', 'locus_index', 'n_het', 'r', 'p'])
tf_df.to_csv(out_tfs_loci, sep='\t', index=False)



# ''' save tf and ase expression in a csv file '''
# if(len(tf_loc_pairs)>0):
#     li = tf_loc_pairs[0][1]
#     tf_ase = [(e,a) for h,e,a 
#                        in zip(het_data.loc[het_data.index, het_data.columns[li]], 
#                               tf_expr_data.loc[het_data.index, tfg],
#                               ase_data.loc[het_data.index,ase_data.columns[li]]) 
#                        if h==1]
#     tf_ase_df = pd.DataFrame(tf_ase, columns=['tf_expr', 'ase_expr'])
#     tf_ase_df.to_csv('results/gene-tf-ase-data.csv')
#     x = tf_ase_df['tf_expr'].values
#     y = tf_ase_df['ase_expr'].values
#     plt.scatter(x, y, alpha=0.5)
