''' 
This script creates histograms of ase-covariate correlations. 
Here, only significantly correlated ase-cov pairs found after multiple-tests correction are used.
'''

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import itertools


''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_covariate_sig_corr_fdr_bh.txt'

# output file
dest_fig_path = 'results/ase_cov_fdr_sig_corr_hist.png'
dest_frequent_cov_path = 'results/ase_cov_fdr_sig_corr_frequent_cov.txt'
# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
corr_data = pd.read_table(cov_ase_corr_data_path, sep='\t', header=0, index_col=False)

# generate correlation histogram
r_vals = corr_data['r'].tolist()
r_vals = [r for r in r_vals if not np.isnan(r)]

# graph settings
nrow = 1
ncol = 1
nbins = 20
xlimit = (-0.75,0.75)
xticks = np.arange(xlimit[0], xlimit[1]+0.1, 0.25)


plt.subplot(nrow, ncol, 1)
h = plt.hist(r_vals, bins=nbins)
plt.axis('auto')
#plt.xlim(xlimit)
#plt.xticks(xticks)

plt.suptitle('Histogram of significant covariate-asesite r')
plt.savefig(dest_fig_path)
plt.close()


covs = corr_data['covariate'].tolist()	
covs.sort()
cov_freq = [(key, len(list(group))) for key, group in itertools.groupby(covs)]
cov_freq.sort(key=lambda x: x[1], reverse=True);
with open(dest_frequent_cov_path, 'w') as outFile:
	outFile.write('\n'.join([str(cvf[0]) + '\t' + str(cvf[1]) for cvf in cov_freq]))
