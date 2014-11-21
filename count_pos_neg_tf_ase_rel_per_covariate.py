''' 
This script counts number of positive and negative tf-ase relations for every significant covariates.
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
frequent_cov_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_cov_fdr_sig_corr_frequent_cov.txt'
tf_ase_data_path = 'results/TF-ASE-Correlation-2014-09-16/significant_tf_ase_fdr_bh.txt'

# output file
dest_path = 'results/cov_pos_neg_r_count.txt'

# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
bio_corr_data = pd.read_table(bio_cov_sig_ase_corr_data_path, sep='\t', header=0, index_col=False)
tech_corr_data = pd.read_table(tech_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)
frequent_cov_data = pd.read_table(frequent_cov_data_path, sep='\t', header=None, index_col=0)
tf_ase_data = pd.read_table(tf_ase_data_path, sep='\t', header=0, index_col=False)

# merge biological and technical data
corr_data = pd.concat([bio_corr_data, tech_corr_data])

print('counting positive and negative tf-ase correlation per covariate')
tf_ase_data['tf_ase_pos'] = tf_ase_data['r'] >= 0
groups = tf_ase_data.groupby(['locus_index', 'tf_ase_pos']).groups
groupKeys = groups.keys()

with open(dest_path, 'w') as outFile:
    for i in range(frequent_cov_data.shape[0]):
        fc = frequent_cov_data.index.values[i]
        posCount = len(groups[(fc,True)]) if (fc, True) in groupKeys else 0
        negCount = len(groups[(fc,False)]) if (fc, False) in groupKeys else 0
        outFile.write('\t'.join([str(fc), str(posCount), str(negCount)]) + '\n');























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

p.savefig(dest_fig, dpi=300)
p.close()

