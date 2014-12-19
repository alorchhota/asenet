import os
import pandas as pd

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
ase_loc_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ase_locus_annot.txt'
bias_info_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/processed_data/snparray.autosomal.uniq.bias.txt'

bias_dest_path = 'results/bias.txt'

# set working directory
os.chdir(home_dir)

print('reading dgn ase locus data ...')
ase_loc_ann_data = pd.read_table(ase_loc_path, sep=' ', header=0, index_col=None)
#ase_loc_ann_data = ase_loc_ann_data.iloc[:100, :]

print('reading dgn bias data ...')
bias_data = pd.read_table(bias_info_path, sep='\t', header=0, index_col=None)
bias_cols = ['chr', 'cord', 'Percentage_ref']
bias_data = bias_data.loc[:, bias_cols]

merged_data = pd.merge(ase_loc_ann_data, bias_data, left_on=['chr', 'asePos'],  right_on=['chr', 'cord'], sort=False, how='inner')

# check if every ase_locus is present in merged data
if merged_data.shape[0] != ase_loc_ann_data.shape[0]:
    print('some ase loci do not have bias information.')
else:
    print('creating bias information ...')
    # hashing first to ensure data in the exact same order as ase loci
    ase_bias_dict = {str(int(row[1]['chr'])) + ':' + str(int(row[1]['asePos'])) : int(row[1]['Percentage_ref'] == 0.5) for row in merged_data.iterrows()}
    ase_bias = [str(int(ase_bias_dict[str(int(row[1]['chr'])) + ':' + str(int(row[1]['asePos']))])) for row in ase_loc_ann_data.iterrows()]

    print('saving bias information ...')
    text = '\n'.join(ase_bias)
    with open(bias_dest_path, 'w') as outFile:
        outFile.write(text)

    print('Done. See output in File: ' + bias_dest_path)