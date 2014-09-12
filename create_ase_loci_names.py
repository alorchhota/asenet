import os
import pandas as pd

''' setting variables '''
home_dir = '/home/asaha6/github/asenet'
gencode_annot_path = 'data/gencode.v19.annotation.processed.gtf'
ase_loc_path = '/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ase_locus_annot.txt'


ase_loc_name_dest_path = 'results/ase_locus_name.txt'

# set working directory
os.chdir(home_dir)

# import files from current directory
import gtf

print('reading gencode annotations ...')
gen_annot = gtf.GencodeGTFReader(gencode_annot_path, processed=True)
gen_annot.data = gen_annot.data.loc[:, ['chr', 'feature', 'start', 'end', 'gene_name']]

print('reading dgn ase locus data ...')
ase_loc_ann_data = pd.read_table(ase_loc_path, sep=' ', header=0, index_col=None)

# reduce the size
#ase_loc_ann_data = ase_loc_ann_data.iloc[0:20,:] 

print('extracting locus gene names ...')
locus_genes = [gen_annot.get_gene_name_from_locus(locus[0], locus[1]) for index, locus in ase_loc_ann_data.iterrows()]
 
print('saving locus gene names ...')
with open(ase_loc_name_dest_path, 'w') as outFile:
    outFile.write('\n'.join([ str(g) for g in locus_genes]))
