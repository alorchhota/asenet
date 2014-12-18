import os
import pandas as pd
import numpy as np

''' setting variables chaning something...'''
home_dir = '/home/asaha6/github/asenet'

# input
sig_tf_ase_path = 'results/TF-ASE-Correlation-2014-12-17/significant_tf_ase_fdr_bh.txt'
sig_tf_tf_path = 'results/TF-TF-Correlations-2014-09-18/tf_tf_test_statistics_2014-09-18_11-28-38_fdr_bh.txt'

# output
triplet_dest_path = 'results/sig_tf_tf_ase_both_fdr.txt'


# set working directory
os.chdir(home_dir)

print('reading test statistics ...')
tf_ase_data = pd.read_table(sig_tf_ase_path, sep='\t', header=0, index_col=None)
tf_tf_data = pd.read_table(sig_tf_tf_path, sep='\t', header=0, index_col=None)


print('creating a map from a tf to ASEs')
groups = tf_ase_data.groupby(['tf']).groups
tf_ase_map = {}
for tf in groups:
	ase_df = tf_ase_data.iloc[groups[tf],]
	ases = list(ase_df.loc[:,'gene'])
	tf_ase_map[tf] = ases

#print(tf_ase_map)

print('finding tf-tf-ase triplets...')
# include if any TF is associated with an ASE
'''
fh = open(triplet_dest_path, "w+")
fh.write('\t'.join(['tf1','tf2','gene', 'tf', 'tf_tf_r', 'tf_tf_p', 'tf_gene_r', 'tf_gene_p']) + '\n');
for tf_tf in tf_tf_data.iterrows():
    tf1 = tf_tf[1]['tf1']
    tf2 = tf_tf[1]['tf2']
    tf_tf_r = tf_tf[1]['r']
    tf_tf_p = tf_tf[1]['p']
    #
    if tf1 in groups:
    	ase_df = tf_ase_data.iloc[groups[tf1],]
    	for tf_ase in ase_df.iterrows():
    		tf = tf_ase[1]['tf']
    		gene = tf_ase[1]['gene']
    		tf_gene_r = tf_ase[1]['r']
    		tf_gene_p = tf_ase[1]['p']
    		fh.write('\t'.join([tf1, tf2, gene, tf, str(tf_tf_r), str(tf_tf_p), str(tf_gene_r), str(tf_gene_p)]) + '\n')
	#
    if tf2 in groups:
    	ase_df = tf_ase_data.iloc[groups[tf2],]
    	for tf_ase in ase_df.iterrows():
    		tf = tf_ase[1]['tf']
    		gene = tf_ase[1]['gene']
    		tf_gene_r = tf_ase[1]['r']
    		tf_gene_p = tf_ase[1]['p']
    		fh.write('\t'.join([tf1, tf2, gene, tf, str(tf_tf_r), str(tf_tf_p), str(tf_gene_r), str(tf_gene_p)]) + '\n')


fh.close();
'''


# include if both TFs are associated with same ASE
fh = open(triplet_dest_path, "w+")
fh.write('\t'.join(['tf1','tf2','gene', 'tf1_tf2_r', 'tf1_gene_r', 'tf2_gene_r']) + '\n');
for tf_tf in tf_tf_data.iterrows():
    tf1 = tf_tf[1]['tf1']
    tf2 = tf_tf[1]['tf2']
    tf_tf_r = tf_tf[1]['r']
    tf_tf_p = tf_tf[1]['p']
    #
    if tf1 in groups and tf2 in groups:
    	common_genes = set(tf_ase_map[tf1]) & set(tf_ase_map[tf2])
    	ase1_df = tf_ase_data.iloc[groups[tf1],]
    	ase2_df = tf_ase_data.iloc[groups[tf2],]
    	for g in common_genes:
    		# find rows in ase1_df with g
    		tf1_gene_r = []
    		for row in ase1_df.iterrows():
    			if row[1]['gene'] == g:
    				tf1_gene_r.append(row[1]['r'])
    		# find rows in ase2_df with g
    		tf2_gene_r = []
    		for row in ase2_df.iterrows():
    			if row[1]['gene'] == g:
    				tf2_gene_r.append(row[1]['r'])
    		# report all combinations of rows1 and rows2 
    		if len(tf1_gene_r) > 0 and len(tf2_gene_r) > 0:
    			comb = [(tf1, tf2, g, tf_tf_r, r1, r2) for r1 in tf1_gene_r for r2 in tf2_gene_r]
    			text = '\n'.join(['\t'.join([str(item) for item in c]) for c in comb]) + '\n'
    			fh.write(text)
    
    
fh.close()
