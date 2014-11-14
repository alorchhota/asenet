''' 
This script creates histograms of ase-covariate correlations. 
Here, NOT all ASEs, but only ASEs appeared in significant TF-ASE correlation pairs have been used.
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
#cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/sig_ase_biological_covariate_corr.txt'
cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/sig_ase_technical_covariate_corr.txt'

# output file
#dest_fig_path = 'results/cov_bio_sig_ase_corr_hist.png'
dest_fig_path = 'results/cov_tech_sig_ase_corr_hist.png'

# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
corr_data = pd.read_table(cov_ase_corr_data_path, sep='\t', header=0, index_col=False)

# generate correlation histogram
r_vals = corr_data['r'].tolist()
r_vals = [r for r in r_vals if not np.isnan(r)]
r_vals_tf_ase_pos = [ row['r']  for idx, row in corr_data.iterrows() if row['tf_ase_r']>=0 and not np.isnan(row['r'])]
r_vals_tf_ase_neg = [ row['r']  for idx, row in corr_data.iterrows() if row['tf_ase_r']<0 and not np.isnan(row['r'])]

#ks_stat = stats.ks_2samp(r_vals_tf_ase_pos, r_vals_tf_ase_neg)
#print(ks_stat)

xlimit = (-0.75,0.75)
xticks = np.arange(xlimit[0], xlimit[1]+0.1, 0.25)
ylimit = 300
nbins = 20

plt.subplot(2, 2, 1)
h = plt.hist(r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.xticks(xticks)
plt.title('All r')
#plt.xlabel('r')

plt.subplot(2, 2, 3)
h = plt.hist(r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.ylim((0,100))
plt.xticks(xticks)
plt.title('All r (zoomed)')
#plt.xlabel('r')

plt.subplot(2, 2, 2)
h = plt.hist(r_vals_tf_ase_pos, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.ylim((0,ylimit))
plt.xticks(xticks)
plt.title('For positive tf-ase corr')
#plt.xlabel('r')

plt.subplot(2, 2, 4)
h = plt.hist(r_vals_tf_ase_neg, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.ylim((0,ylimit))
plt.xticks(xticks)
plt.title('For negative tf-ase corr')
#plt.xlabel('r')

plt.suptitle('Histogram of covariate-asesite r')
plt.savefig(dest_fig_path)
plt.close()
