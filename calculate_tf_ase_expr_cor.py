import pandas as pd
from scipy import stats
import numpy as np


tfAsePairFilepath = 'results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_fdr_bh.txt'
validity_data_path = 'data/ase_validity.txt'
expr_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/data_used_for_eqtl_study/trans_data.txt'

dest_path ='results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_fdr_bh_exprCor.txt'


print('reading validity data ...')
validity_data = pd.read_table(validity_data_path, sep='\t', header=None, index_col=0)

print('reading expresson data ...')
expr_data = pd.read_table(expr_data_path, sep='\t', header=0, index_col=0)
expr_data = expr_data.loc[validity_data.index.values, :]
expr_genes = list(expr_data.columns)

print('reading tf-ase pair data ...')
tf_ase_pair_data = pd.read_table(tfAsePairFilepath, sep='\t', header=0)

print('calculating tf expresson - ase gene expression correlations ...')
tf_ase_expr_cor_r_all = []
tf_ase_expr_cor_p_all = []
tf_ase_expr_cor_r_valid = []
tf_ase_expr_cor_p_valid = []

for index, r in tf_ase_pair_data.iterrows():
    tf = r['tf']
    g = r['gene']
    idx = r['locus_index']
    if tf in expr_genes and g in expr_genes:
        tf_expr = expr_data.loc[:,tf]
        g_expr = expr_data.loc[:,g]
        # calculate correlation using all samples
        cor = stats.spearmanr(tf_expr, g_expr)
        tf_ase_expr_cor_r_all.append(cor[0])
        tf_ase_expr_cor_p_all.append(cor[1])
        # calculate correlation using valid samples only
        validIndexes = np.where(validity_data.iloc[:,idx]==1)[0]
        tf_expr = tf_expr.iloc[validIndexes]
        g_expr = g_expr.iloc[validIndexes]
        cor = stats.spearmanr(tf_expr, g_expr)
        tf_ase_expr_cor_r_valid.append(cor[0])
        tf_ase_expr_cor_p_valid.append(cor[1])
        #print(r['tf'] + ' ' + r['gene'] + ' ' + str(cor[0]))
    else:
        tf_ase_expr_cor_r_all.append('-')
        tf_ase_expr_cor_p_all.append('-')
        tf_ase_expr_cor_r_valid.append('-')
        tf_ase_expr_cor_p_valid.append('-')
    

tf_ase_pair_data['tf_ase_expr_cor_r_all'] = tf_ase_expr_cor_r_all
tf_ase_pair_data['tf_ase_expr_cor_p_all'] = tf_ase_expr_cor_p_all
tf_ase_pair_data['tf_ase_expr_cor_r_valid'] = tf_ase_expr_cor_r_valid
tf_ase_pair_data['tf_ase_expr_cor_p_valid'] = tf_ase_expr_cor_p_valid

tf_ase_pair_data.to_csv(path_or_buf=dest_path, sep='\t', header=True, index=False)

