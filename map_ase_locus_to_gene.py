import os
import pandas as pd
import argparse
import numpy as np

''' argument parsing '''
print('parsing arguments ...')
parser = argparse.ArgumentParser()
parser.add_argument('-home',
                    help='home directory',
                    default='.', )
parser.add_argument('-genannotdata',
                    help='path to gencode annotation data (processed) file.',
                    default='data/gencode.v19.annotation.processed.gtf')
parser.add_argument('-lociannotdata',
                    help='path to ase loci annotation data file.',
                    default='/scratch1/langmead-fs1/data/big_public_datasets/dgn/ase/ase_locus_annot.txt')
parser.add_argument('-annotdest',
                    help='output path to annotated tf-asesite pairs.',
                    default='results/ase_locus_gene_annot.txt')
args = parser.parse_args()

''' setting variables '''
home_dir = args.home
gencode_annot_path = args.genannotdata
ase_loc_path = args.lociannotdata
annot_dest_path = args.annotdest


# set working directory
os.chdir(home_dir)

# import files from current directory

import gtf

print('reading gencode annotations ...')
gen_annot = gtf.GencodeGTFReader(gencode_annot_path, processed=True)
gen_annot.data = gen_annot.data.loc[:, ['chr', 'feature', 'start', 'end', 'gene_name', 'gene_type']]

print('reading dgn ase locus data ...')
ase_loc_ann_data = pd.read_table(ase_loc_path, sep=' ', header=0, index_col=None)
n_loci = ase_loc_ann_data.shape[0]


def isPseudogene(gene_type):
    return(gene_type.endswith('pseudogene'))

def containsPseudogene(gene_types):
    return(any([isPseudogene(gt) for gt in gene_types]))

# list to store annotations of all loci
annotation = ['' for _ in range(n_loci)]
pseudoStatus = [False for _ in range(n_loci)]
n_chr = 23

for chr in range(1, n_chr+1):
    print('chr#' + str(chr))
    # find gencode annotation data for chromosome chr
    matched_chr = gen_annot.data.loc[:, 'chr'] == ('chr'+str(chr))
    matched_index = np.where(matched_chr)[0]
    chr_annotdata = gen_annot.data.iloc[matched_index, :]
    #
    # find locus positions in chromosome chr
    chr_locusdata = [(idx, row['chr'], row['asePos']) for idx, row in ase_loc_ann_data.iterrows() if row['chr']==chr]
    #
    # annotate all loci in chr
    for ld in chr_locusdata:
        pos = ld[2]
        #
        # find genes at pos
        matched_start = chr_annotdata.loc[:, 'start'] <= pos
        matched_end = chr_annotdata.loc[:, 'end'] >= pos
        matched = matched_start & matched_end
        matched_index = np.where(matched)[0]
        matched_entities = chr_annotdata.iloc[matched_index, :]
        gene_names = set(matched_entities.loc[:, 'gene_name'])
        gene_types = set(matched_entities.loc[:, 'gene_type'])
        #
        # if multiple genes mapped, take the expressed one
        gene_name = ''
        ps = False
        if len(gene_names) == 0:
            gene_name = '-'
        elif len(gene_names) == 1:
            gene_name = gene_names.pop()
            ps = isPseudogene(matched_entities.iloc[0, 5])
            if len(gene_types) > 1:
                print('Multiple gene types: ' + str(ld[0]))
        else:
            is_exon = matched_entities.loc[:, 'feature'] == 'exon'
            if sum(is_exon) > 0:
                exon_index = np.where(is_exon)[0]
                gene_name = matched_entities.iloc[exon_index[0], 4] # take first gene_name
                ps = isPseudogene(matched_entities.iloc[exon_index[0], 5])
            else:
                gene_name = matched_entities.iloc[0, 4] # take first gene_name
                ps = isPseudogene(matched_entities.iloc[0, 5])
        #
        # store annotations
        annotation[ld[0]] = gene_name
        pseudoStatus[ld[0]] = ps

print('saving annotations into file:' + annot_dest_path + ' ...')
with open(annot_dest_path, 'w') as fh:
    text = '\n'.join([pair[0] + '\t' + str(int(pair[1])) for pair in zip(annotation, pseudoStatus)])
    fh.write(text)

print('Done!' )