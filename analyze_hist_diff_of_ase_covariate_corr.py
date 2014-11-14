''' 
This script calcualtes difference among histograms of ase-covariate correlations. 
'''


import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats




''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
bio_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/sig_ase_biological_covariate_corr.txt'
tech_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/sig_ase_technical_covariate_corr.txt'

# output file
dest_path = 'results/cov_ase_hist_diff.txt'

# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
bio_corr_data = pd.read_table(bio_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)
tech_corr_data = pd.read_table(tech_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)


print('analyzing biological covariate - ase correlation histograms...')
bio_r_vals = bio_corr_data['r'].tolist()
bio_r_vals = [r for r in bio_r_vals if not np.isnan(r)]
bio_r_vals_tf_ase_pos = [ row['r']  for idx, row in bio_corr_data.iterrows() if row['tf_ase_r']>=0 and not np.isnan(row['r'])]
bio_r_vals_tf_ase_neg = [ row['r']  for idx, row in bio_corr_data.iterrows() if row['tf_ase_r']<0 and not np.isnan(row['r'])]

with open(dest_path, 'w') as outFile:
    outFile.write( '\t'.join(['cov_type', 'comparision', 'D', 'p']) + '\n')
    # all (background) vs positive tf-ase-r
    bio_ks_stat = stats.ks_2samp(bio_r_vals, bio_r_vals_tf_ase_pos)
    outFile.write('\t'.join(['biological', 'bg-vs-pos', str(bio_ks_stat[0]), str(bio_ks_stat[1])]) + '\n')
    # all (background) vs negative tf-ase-r
    bio_ks_stat = stats.ks_2samp(bio_r_vals, bio_r_vals_tf_ase_neg)
    outFile.write('\t'.join(['biological', 'bg-vs-neg', str(bio_ks_stat[0]), str(bio_ks_stat[1])]) + '\n')
    # positive vs. negative
    bio_ks_stat = stats.ks_2samp(bio_r_vals_tf_ase_pos, bio_r_vals_tf_ase_neg)
    outFile.write('\t'.join(['biological', 'pos-vs-neg', str(bio_ks_stat[0]), str(bio_ks_stat[1])]) + '\n')
    

print('analyzing biological covariate - ase correlation histograms...')
tech_r_vals = tech_corr_data['r'].tolist()
tech_r_vals = [r for r in tech_r_vals if not np.isnan(r)]
tech_r_vals_tf_ase_pos = [ row['r']  for idx, row in tech_corr_data.iterrows() if row['tf_ase_r']>=0 and not np.isnan(row['r'])]
tech_r_vals_tf_ase_neg = [ row['r']  for idx, row in tech_corr_data.iterrows() if row['tf_ase_r']<0 and not np.isnan(row['r'])]

with open(dest_path, 'a') as outFile:
    # all (background) vs positive tf-ase-r
    tech_ks_stat = stats.ks_2samp(tech_r_vals, tech_r_vals_tf_ase_pos)
    outFile.write('\t'.join(['technical', 'bg-vs-pos', str(tech_ks_stat[0]), str(tech_ks_stat[1])]) + '\n')
    # all (background) vs negative tf-ase-r
    tech_ks_stat = stats.ks_2samp(tech_r_vals, tech_r_vals_tf_ase_neg)
    outFile.write('\t'.join(['technical', 'bg-vs-neg', str(tech_ks_stat[0]), str(tech_ks_stat[1])]) + '\n')
    # positive vs. negative
    tech_ks_stat = stats.ks_2samp(tech_r_vals_tf_ase_pos, tech_r_vals_tf_ase_neg)
    outFile.write('\t'.join(['technical', 'pos-vs-neg', str(tech_ks_stat[0]), str(tech_ks_stat[1])]) + '\n')
    
    
print('analyzing biological+technical covariate - ase correlation histograms...')
r_vals = bio_r_vals + tech_r_vals
r_vals_tf_ase_pos = bio_r_vals_tf_ase_pos + tech_r_vals_tf_ase_pos
r_vals_tf_ase_neg = bio_r_vals_tf_ase_neg + tech_r_vals_tf_ase_neg

with open(dest_path, 'a') as outFile:
    # all (background) vs positive tf-ase-r
    ks_stat = stats.ks_2samp(r_vals, r_vals_tf_ase_pos)
    outFile.write('\t'.join(['bio+tech', 'bg-vs-pos', str(ks_stat[0]), str(ks_stat[1])]) + '\n')
    # all (background) vs negative tf-ase-r
    ks_stat = stats.ks_2samp(r_vals, r_vals_tf_ase_neg)
    outFile.write('\t'.join(['bio+tech', 'bg-vs-neg', str(ks_stat[0]), str(ks_stat[1])]) + '\n')
    # positive vs. negative
    ks_stat = stats.ks_2samp(r_vals_tf_ase_pos, r_vals_tf_ase_neg)
    outFile.write('\t'.join(['bio+tech', 'pos-vs-neg', str(ks_stat[0]), str(ks_stat[1])]) + '\n')
    
