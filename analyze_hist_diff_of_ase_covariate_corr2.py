''' 
This script visualizes different histograms of ase-covariate correlations to check how much they are different from background.
Here, background means all possible ase-covariate correlation(r).
Two histograms are compared with background. 
1) ase-cov correlation where tf-ase-r is positive (non-negative)
2) ase-cov correlation where tf-ase-r is negative
'''


import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as p
from scipy import stats




''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
bio_cov_sig_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/sig_ase_biological_covariate_corr.txt'
tech_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/sig_ase_technical_covariate_corr.txt'
bio_cov_all_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_biological_covariate_corr.txt'
tech_cov_all_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_technical_covariate_corr.txt'


# output file
dest_fig = 'results/cov_ase_hist_line.png'

# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
bio_corr_data = pd.read_table(bio_cov_sig_ase_corr_data_path, sep='\t', header=0, index_col=False)
tech_corr_data = pd.read_table(tech_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)

bio_corr_data_all = pd.read_table(bio_cov_all_ase_corr_data_path, sep='\t', header=0, index_col=False)
tech_corr_data_all = pd.read_table(tech_cov_all_ase_corr_data_path, sep='\t', header=0, index_col=False)

## graph configurations
nbins = 20
nrow = 2
ncol = 2

print('analyzing biological covariate - ase correlation histograms...')
bio_r_vals = bio_corr_data['r'].tolist()
bio_r_vals = [r for r in bio_r_vals if not np.isnan(r)]
bio_r_vals_tf_ase_pos = [ row['r']  for idx, row in bio_corr_data.iterrows() if row['tf_ase_r']>=0 and not np.isnan(row['r'])]
bio_r_vals_tf_ase_neg = [ row['r']  for idx, row in bio_corr_data.iterrows() if row['tf_ase_r']<0 and not np.isnan(row['r'])]

bio_r_vals_all = bio_corr_data_all['r'].tolist()
bio_r_vals_all = [r for r in bio_r_vals_all if not np.isnan(r)]


p.subplot(nrow, ncol, 1)
y,binEdges=np.histogram(bio_r_vals_all, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='bg', color='blue')

#y,binEdges=np.histogram(bio_r_vals, bins=nbins, normed=True)
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#p.plot(bincenters,y,'-', label='tf-ase')

y,binEdges=np.histogram(bio_r_vals_tf_ase_pos, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(+r)', color='green')

y,binEdges=np.histogram(bio_r_vals_tf_ase_neg, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(-r)', color='red')

#p.legend()
p.title('Biological Covariates')
p.show()


print('analyzing technical covariate - ase correlation histograms...')
tech_r_vals = tech_corr_data['r'].tolist()
tech_r_vals = [r for r in tech_r_vals if not np.isnan(r)]
tech_r_vals_tf_ase_pos = [ row['r']  for idx, row in tech_corr_data.iterrows() if row['tf_ase_r']>=0 and not np.isnan(row['r'])]
tech_r_vals_tf_ase_neg = [ row['r']  for idx, row in tech_corr_data.iterrows() if row['tf_ase_r']<0 and not np.isnan(row['r'])]

tech_r_vals_all = tech_corr_data_all['r'].tolist()
tech_r_vals_all = [r for r in tech_r_vals_all if not np.isnan(r)]


p.subplot(nrow, ncol, 2)
y,binEdges=np.histogram(tech_r_vals_all, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='bg', color='blue')

#y,binEdges=np.histogram(tech_r_vals, bins=nbins, normed=True)
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#p.plot(bincenters,y,'-', label='tf-ase')

y,binEdges=np.histogram(tech_r_vals_tf_ase_pos, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(+r)', color='green')

y,binEdges=np.histogram(tech_r_vals_tf_ase_neg, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(-r)', color='red')

#p.legend()
p.title('Technical Covariates')
p.show()



print('analyzing biological+technical covariate - ase correlation histograms...')
r_vals = bio_r_vals + tech_r_vals
r_vals_tf_ase_pos = bio_r_vals_tf_ase_pos + tech_r_vals_tf_ase_pos
r_vals_tf_ase_neg = bio_r_vals_tf_ase_neg + tech_r_vals_tf_ase_neg
r_vals_all = bio_r_vals_all + tech_r_vals_all 

p.subplot(nrow, ncol, 3)
y,binEdges=np.histogram(r_vals_all, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='bg', color='blue')

#y,binEdges=np.histogram(r_vals, bins=nbins, normed=True)
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#p.plot(bincenters,y,'-', label='tf-ase')

y,binEdges=np.histogram(r_vals_tf_ase_pos, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(+r)', color='green')

y,binEdges=np.histogram(r_vals_tf_ase_neg, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(-r)', color='red')

#p.legend()
p.title('Both Covariates')
p.show()


p.subplot(nrow, ncol, 4)
dummy = [0,0,0,0,0]
y,binEdges=np.histogram(dummy, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='bg', color='blue')

#y,binEdges=np.histogram(r_vals, bins=nbins, normed=True)
#bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#p.plot(bincenters,y,'-', label='tf-ase')

y,binEdges=np.histogram(dummy, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(+r)', color='green')

y,binEdges=np.histogram(dummy, bins=nbins, normed=True)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
p.plot(bincenters,y,'-', label='tf-ase(-r)', color='red')

p.legend()
p.title('Legends')
p.show()

p.savefig(dest_fig, dpi=300)
p.close()

