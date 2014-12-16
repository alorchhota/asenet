import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

''' setting variables chaning something...'''
home_dir = '/home/asaha6/github/asenet'

# input
sig_tf_ase_path_noabs = 'results/TF-ASE-Correlation-2014-12-15/significant_tf_ase_fdr_bh.txt'
sig_tf_ase_path_abs = 'results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_fdr_bh.txt'
tf_expr_data_path = 'data/tf_expr_data.txt'
validity_data_path_noabs = 'data/ase_validity_noabs.txt'
validity_data_path_abs = 'data/ase_validity.txt'
ase_data_path_noabs = 'data/ase_noabs.txt'
ase_data_path_abs = 'data/ase.txt'
ref_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ref.txt'
alt_data_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/alt.txt'

zero_cross_thresh = 0.1 # 10%

# output
hit_miss_dest_path = 'results/hit_missed_without_abs.txt'
fig_dest_dir = 'results/TF-ASE-Correlation-2014-12-15/fig_hit_missed_without_abs/'

# set working directory
os.chdir(home_dir)

print('reading data ...')
tf_ase_df_abs = pd.read_table(sig_tf_ase_path_abs, sep='\t', header=0, index_col=None)
tf_ase_df_noabs = pd.read_table(sig_tf_ase_path_noabs, sep='\t', header=0, index_col=None)

ase_df_abs = pd.read_table(ase_data_path_abs, sep='\t', header=None, index_col=0)
ase_df_noabs = pd.read_table(ase_data_path_noabs, sep='\t', header=None, index_col=0)

validity_df_abs = pd.read_table(validity_data_path_abs, sep='\t', header=None, index_col=0)
validity_df_noabs = pd.read_table(validity_data_path_noabs, sep='\t', header=None, index_col=0)

ref_data = pd.read_table(ref_data_path, sep=' ', header=None, index_col=0)
alt_data = pd.read_table(alt_data_path, sep=' ', header=None, index_col=0)

tf_expr_data = pd.read_table(tf_expr_data_path, sep='\t', header=0, index_col=0)



print('comparing tf-ase pairs with and without absolute values ...')
tf_ase_pairs_abs = [(ta[1]['tf'], ta[1]['locus_index'], ta[1]['gene']) for ta in tf_ase_df_abs.iterrows()]
tf_ase_pairs_noabs = [(ta[1]['tf'], ta[1]['locus_index'], ta[1]['gene']) for ta in tf_ase_df_noabs.iterrows()]

abs_pairs = set(tf_ase_pairs_abs)
noabs_pairs = set(tf_ase_pairs_noabs)
common_pairs = abs_pairs.intersection(noabs_pairs)
missed_pairs = abs_pairs.difference(noabs_pairs)
new_pairs = noabs_pairs.difference(abs_pairs)

with open(hit_miss_dest_path, 'w') as outFile:
    outFile.write('#common_hits:\t' + str(len(common_pairs)) + '\n')
    outFile.write('#missed_hits:\t' + str(len(missed_pairs)) + '\n')
    outFile.write('#new_hits:\t' + str(len(new_pairs)) + '\n')
    outFile.write('\nMissed Hits:\n')
    outFile.write('tf\tlocus_index\tgene\n')
    text = '\n'.join(['\t'.join([str(item) for item in mis]) 
                      for mis in missed_pairs])
    outFile.write(text)


print('plotting missed hits ...')
def drawPlot(tf, l, l_gene):
    v_index = np.where(validity_df_abs.iloc[:,l]==1)[0]
    tf_expr = tf_expr_data.iloc[v_index, :]
    tf_expr = tf_expr.loc[:,tf]
    ase_expr_abs = ase_df_abs.iloc[v_index, l]
    ase_expr_noabs = ase_df_noabs.iloc[v_index, l]
    ref_count = ref_data.iloc[v_index, l]
    alt_count = alt_data.iloc[v_index,l]
    #print('plot...')
    f, subplots = plt.subplots(2,2)
    # |ase| vs tf
    gi,gj=0,0
    sp = subplots[gi,gj].scatter(tf_expr, ase_expr_abs)
    cor = stats.spearmanr(tf_expr, ase_expr_abs)
    subplots[gi,gj].legend([sp],['|ase|'])
    subplots[gi,gj].set_ylabel('Locus: ' + l_gene + '(' + str(l)+')')
    subplots[gi,gj].set_title('r:%.2f, p:%.2e'%(cor[0],cor[1]))
    #subplots[gi,gj].set_title('ASE vs TF Expression')
    # ase vs tf
    gi,gj=0,1
    sp = subplots[gi,gj].scatter(tf_expr, ase_expr_noabs)
    subplots[gi,gj].legend([sp],['ase'])
    cor = stats.spearmanr(tf_expr, ase_expr_noabs)
    subplots[gi,gj].set_title('r:%.2f, p:%.2e'%(cor[0],cor[1]))
    # #ref vs tf
    gi,gj = 1,0
    sp = subplots[gi,gj].scatter(tf_expr, ref_count)
    subplots[gi,gj].legend([sp],['#ref_reads'])
    subplots[gi,gj].set_xlabel('TF: ' + tf)
    # #alt vs tf
    gi,gj = 1,1
    sp = subplots[gi,gj].scatter(tf_expr, alt_count)
    subplots[gi,gj].legend([sp],['#alt_reads'])
    plt.savefig(fig_dest_dir + "/" + tf + "_" + str(l) + ".png")
    plt.close()
    return
    
#drawPlot(tf, l, l_gene, r_abs, p_abs, r_noabs, p_noabs):
for mis in missed_pairs:
    drawPlot(mis[0],mis[1],mis[2])
    
print('finding zero crossings ...')
zero_crossings = []
for tap in noabs_pairs:
    l = tap[1]
    v_index = np.where(validity_df_noabs.iloc[:,l]==1)[0]
    ase_expr = ase_df_noabs.iloc[v_index, l]
    pos_count = sum([1 for e in ase_expr if e>0])
    neg_count = len(ase_expr) - pos_count
    zero_cross = min([pos_count, neg_count])/(len(ase_expr)+0.0)
    if zero_cross > zero_cross_thresh:
        zero_crossings.append(tap)

with open(hit_miss_dest_path, 'a') as outFile:
    outFile.write('\n\n#zero_crossing:\t' + str(len(zero_crossings)) + '\n')
    outFile.write('\n\nZero Crossings:\n')
    outFile.write('tf\tlocus_index\tgene\n')
    text = '\n'.join(['\t'.join([str(item) for item in zc]) 
                      for zc in zero_crossings])
    outFile.write(text)
    
