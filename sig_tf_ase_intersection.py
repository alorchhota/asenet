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
parser.add_argument('-data1',
                    help='path to significant TF-ASE pairs data file 1.',
                    default='/home/ashis/temp/asenet_results/TF-ASE-Correlation-2015-01-07-pseudo50/significant_tf_ase_bonferroni.txt')
parser.add_argument('-data2',
                    help='path to significant TF-ASE pairs data file 2.',
                    default='/home/ashis/temp/asenet_results/TF-ASE-Correlation-2015-01-07-pseudo20/significant_tf_ase_bonferroni.txt')
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
args = parser.parse_args()


''' setting variables '''
home_dir = args.home
data1_path = args.data1
data2_path = args.data2

# set working directory
os.chdir(home_dir)

print('reading data ...')
data1 = pd.read_table(data1_path, sep='\t', header=0, index_col=None)
data2 = pd.read_table(data2_path, sep='\t', header=0, index_col=None)

# function to generate intersection summary
def intersectSummary(x, y, by):
    if by is None:
        by = set(x.columns.values).intersection(set(y.columns.values))
    #
    x1 =  [tuple(row[col] for col in by) for idx, row in x.iterrows()]
    y1 =  [tuple(row[col] for col in by) for idx, row in y.iterrows()]
    x1 = set(x1)
    y1 = set(y1)
    xy = x1.intersection(y1)
    to_ret = {'xy':len(xy), 'x':len(x), 'y':len(y)}
    return(to_ret)

intsum = intersectSummary(data1, data2, by=['tf','locus_index'])
print(intsum)