''' 
This script creates histograms of ase-covariate correlations. 
Here, every possible ASE-COV pairs are taken care of.
Note: there are two categories of covariates - biological and technical.
'''

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter



''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
bio_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_biological_covariate_corr.txt'
tech_cov_ase_corr_data_path = 'results/COV-ASE-Correlations-2014-11-05/ase_technical_covariate_corr.txt'

# output file
dest_fig_path = 'results/ase_cov_corr_hist.png'

# set working directory
os.chdir(home_dir)

print('reading cov-ase correlation data ...')
bio_corr_data = pd.read_table(bio_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)
tech_corr_data = pd.read_table(tech_cov_ase_corr_data_path, sep='\t', header=0, index_col=False)

# generate correlation histogram
bio_r_vals = bio_corr_data['r'].tolist()
tech_r_vals = tech_corr_data['r'].tolist()
print('#r_vals with nan: ' + str(len(bio_r_vals)) + ", " + str(len(tech_r_vals)))

bio_r_vals = [r for r in bio_r_vals if not np.isnan(r)]
tech_r_vals = [r for r in tech_r_vals if not np.isnan(r)]
print('#r_vals without nan: ' + str(len(bio_r_vals)) + ", " + str(len(tech_r_vals)))

# graph settings
nrow = 3
ncol = 2
nbins = 50
xlimit = (-0.75,0.75)
xticks = np.arange(xlimit[0], xlimit[1]+0.1, 0.25)
ylimit = 500
zoomedYlimit = (0,5000)

# formatter
def to_kilo(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(int(y/1000))
    # The percent symbol needs escaping in latex
    return s + 'k'




plt.subplot(nrow, ncol, 1)
h = plt.hist(bio_r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.xticks(xticks)
plt.ylabel('Biological')
#plt.tick_params(\
#    axis='x',          # changes apply to the x-axis
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    labelbottom='off') # labels along the bottom edge are off
formatter = FuncFormatter(to_kilo)
plt.gca().yaxis.set_major_formatter(formatter)


plt.subplot(nrow, ncol, 2)
h = plt.hist(bio_r_vals, bins=nbins)
plt.xlim(xlimit)
plt.ylim(zoomedYlimit)
plt.xticks(xticks)
#plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
formatter = FuncFormatter(to_kilo)
plt.gca().yaxis.set_major_formatter(formatter)


plt.subplot(nrow, ncol, 3)
h = plt.hist(tech_r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.xticks(xticks)
plt.ylabel('Technical')
#plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
formatter = FuncFormatter(to_kilo)
plt.gca().yaxis.set_major_formatter(formatter)

plt.subplot(nrow, ncol, 4)
h = plt.hist(tech_r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.ylim(zoomedYlimit)
plt.xticks(xticks)
#plt.ylabel('Technical Covariates')
#plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
formatter = FuncFormatter(to_kilo)
plt.gca().yaxis.set_major_formatter(formatter)

plt.subplot(nrow, ncol, 5)
h = plt.hist(bio_r_vals + tech_r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.xticks(xticks)
plt.ylabel('All')
plt.xlabel('r')
formatter = FuncFormatter(to_kilo)
plt.gca().yaxis.set_major_formatter(formatter)

plt.subplot(nrow, ncol, 6)
h = plt.hist(bio_r_vals + tech_r_vals, bins=nbins)
plt.axis('auto')
plt.xlim(xlimit)
plt.ylim(zoomedYlimit)
plt.xticks(xticks)
#plt.title('All Covariates')
plt.xlabel('r')
formatter = FuncFormatter(to_kilo)
plt.gca().yaxis.set_major_formatter(formatter)

plt.suptitle('Histogram of covariate-asesite r')
plt.savefig(dest_fig_path)
plt.close()
