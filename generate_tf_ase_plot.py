import heapq
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

''' argument parsing '''
print('parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
parser.add_argument('-tfdata',
                    help='path to TF expression data file.',
                    default='data/tf_expr_data.txt')
parser.add_argument('-validdata',
                    help='path to validity data file. Only valid loci are included in the analysis.',
                    default='data/ase_validity.txt')
parser.add_argument('-asedata',
                    help='path to ase data file.',
                    default='data/ase.txt')
parser.add_argument('-refdata',
                    help='path to reference read counts data file.',
                    default='/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ref.txt')
parser.add_argument('-altdata',
                    help='path to alternate read counts data file.',
                    default='/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/alt.txt')
parser.add_argument('-plotdir',
                    help='output dir to store plots.',
                    default='results/fig')
parser.add_argument('-sigdata',
                    help='path to annotated significant tf-asesite pairs.',
                    default='results/significant_tf_ase_bonferroni.txt')
args = parser.parse_args()

''' setting variables '''
home_dir = args.home
tf_expr_data_path = args.tfdata
validity_data_path = args.validdata
ase_data_path = args.asedata
ref_data_path = args.refdata
alt_data_path = args.altdata
test_statistics_path = args.sigdata
fig_dest_dir = args.plotdir

# set working directory
os.chdir(home_dir)

print('reading dgn data ...')
ase_data = pd.read_table(ase_data_path, sep='\t', header=None, index_col=0)
validity_data = pd.read_table(validity_data_path, sep='\t', header=None, index_col=0)
ref_data = pd.read_table(ref_data_path, sep=' ', header=None, index_col=0)
alt_data = pd.read_table(alt_data_path, sep=' ', header=None, index_col=0)

print('reading TF expression data ...')
tf_expr_data = pd.read_table(tf_expr_data_path, sep='\t', header=0, index_col=0)

print('reading test statistics ...')
sig_test_stat_data = pd.read_table(test_statistics_path, sep='\t', header=0, index_col=None)



def drawPlot(tf, l, l_gene, r, p):
    v_index = np.where(validity_data.iloc[:,l]==1)[0]
    tf_expr = tf_expr_data.iloc[v_index, :]
    tf_expr = tf_expr.loc[:,tf]
    ase_expr = ase_data.iloc[v_index, l]
    ref_count = ref_data.iloc[v_index, l]
    alt_count = alt_data.iloc[v_index,l]
    #print('plot...')
    f, subplots = plt.subplots(2,2)
    # ase vs tf
    gi,gj=0,0
    sp = subplots[gi,gj].scatter(tf_expr, ase_expr)
    subplots[gi,gj].legend([sp],['ase'])
    subplots[gi,gj].set_ylabel('Locus: ' + l_gene + '(' + str(l)+')')
    #subplots[gi,gj].set_title('ASE vs TF Expression')
    # #samples
    gi,gj=0,1
    subplots[gi,gj].set_xticks([0,1])
    subplots[gi,gj].set_yticks([0,1])
    subplots[gi,gj].text(0.1,0.8, '#samples: ' + str(len(v_index)))
    subplots[gi,gj].text(0.1,0.5, 'r: ' + str(r))
    subplots[gi,gj].text(0.1,0.2, 'p: ' + str(p))
    # #ref vs tf
    gi,gj = 1,0
    sp = subplots[gi,gj].scatter(tf_expr, ref_count)
    subplots[gi,gj].legend([sp],['#ref_reads'])
    subplots[gi,gj].set_xlabel('TF: ' + tf)
    #subplots[gi,gj].set_ylabel('Locus: ' + l_gene + '(' + str(l)+')')
    #subplots[gi,gj].set_title('#RefReads vs TF Expression')
    # #alt vs tf
    gi,gj = 1,1
    sp = subplots[gi,gj].scatter(tf_expr, alt_count)
    subplots[gi,gj].legend([sp],['#alt_reads'])
    #subplots[gi,gj].set_xlabel('TF: ' + tf)
    #subplots[gi,gj].set_ylabel('Locus: ' + l_gene + '(' + str(l)+')')
    #subplots[gi,gj].set_title('#AltReads vs TF Expression')
    plt.savefig(fig_dest_dir + "/" + tf + "_" + str(l) + ".png")
    plt.close()
    return

print('finding top n (highest r) relations')
n = sig_test_stat_data.shape[0]
nthLargestR = heapq.nlargest(n,sig_test_stat_data['r'])[n-1]
sig_idx = np.where(sig_test_stat_data['r']>=nthLargestR)[0]
sig_cor = [(sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'tf'], 
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'locus_index'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'gene'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'r'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'p']) 
           for idx in sig_idx]


print('generating plots ...')
for i in range(n):
    cor = sig_cor[i]
    tf = cor[0]
    l = cor[1]
    l_gene = cor[2]
    r = cor[3]
    p = cor[4]
    print(cor)
    drawPlot(tf, l, l_gene, r, p)

'''
print('finding negative top n (lowest r) relations')
n = 5
nthSmallestR = heapq.nsmallest(n,sig_test_stat_data['r'])[n-1]
sig_idx = np.where(sig_test_stat_data['r']<=nthSmallestR)[0]
sig_cor = [(sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'tf'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'locus_index'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'gene'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'r'],
           sig_test_stat_data.loc[sig_test_stat_data.index[idx], 'p'])
           for idx in sig_idx]


print('generating plots ...')
for i in range(n):
    cor = sig_cor[i]
    tf = cor[0]
    l = cor[1]
    l_gene = cor[2]
    r = cor[3]
    p = cor[4]
    print(cor)
    drawPlot(tf, l, l_gene, r, p)
'''

