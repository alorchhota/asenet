import os
import pandas as pd

''' setting variables '''
home_dir = '/home/ashis/work/github/asenet'
gene_map_path = 'data/gencode.v19.genemap.txt'
het_data_path = 'data/het.txt'
ref_data_path = 'data/ref.txt'
alt_data_path = 'data/alt.txt'
tf_data_path = 'data/tf-nrg2538-s3.txt'
gtex_genes_path = 'data/gtex_genes_used_in_eqtl_study.txt'
gtex_gene_map_by_biomart_org_path = 'data/gtex_genes_map_biomart.org.txt'

# set working directory
os.chdir(home_dir)

# import files from current directory
import tfdb

def read_gene_map(gene_map_path):
    '''create a dictionary from gene map file'''
    cnames = ['gene_id', 'gene_name']
    df = pd.read_table(gene_map_path, sep='\t', header=0, names=cnames)
    id2name_dic = {row['gene_id']:row['gene_name'] for _, row in df.iterrows()}
    return id2name_dic

# generate gene maps'
gene_map = read_gene_map(gene_map_path)
print(gene_map['ENSG00000187223'])

# read TFs
tf = tfdb.TFDB(tf_data_path)
tf_sym = [gene_map[g] if g in gene_map else '-' for g in tf.ensg]

ntf = len(tf.ensg)
combined_tuples = []
for i in range(ntf):
    gid = tf.ensg[i]
    sym_db = tf.sym[i]
    sym_gencode = gene_map[gid] if gid in gene_map else '-'
    combined_tuples.append((gid, sym_db, sym_gencode))

combined_df = pd.DataFrame(combined_tuples)
combined_df.to_csv(path_or_buf='results/map_comparision.txt', sep='\t', header=False)

'''create summary report'''
valid_sym_db = [g for g in tf.sym if g != tfdb.TFDB.na_value]
valid_sym_gencode = [gene_map[gid] for gid in tf.ensg if gid in gene_map]

n_symbols_in_db = len(valid_sym_db)
n_symbols_in_gencode = len(valid_sym_gencode)
n_symbols_agreed_by_both = sum([t[1]==t[2] for t in combined_tuples])

print('#GeneSymbols in PaperDB: ' + str(n_symbols_in_db))
print('#GeneSymbols in Gencode: ' + str(n_symbols_in_gencode))
print('#GeneSymbols agreed by both: ' + str(n_symbols_agreed_by_both))


''' investigate how well these maps do on gtex data '''
gtex_genes = pd.read_csv(gtex_genes_path, names=['gene'])['gene'].tolist()
sym_db_in_gtex = [s for s in valid_sym_db if s in gtex_genes]
sym_gencode_in_gtex = [s for s in valid_sym_gencode if s in gtex_genes]
sym_both_in_gtex = set(sym_db_in_gtex).intersection(set(sym_gencode_in_gtex)) 

print('#GeneSymbols from db mapped do GTEx: ' + str(len(sym_db_in_gtex)))
print('#GeneSymbols from gencode mapped do GTEx: ' + str(len(sym_gencode_in_gtex)))
print('#GeneSymbols from both mapped do GTEx: ' + str(len(sym_both_in_gtex)))


''' check map by biomart.org '''
biomart_org_gene_map_df = pd.read_csv(gtex_gene_map_by_biomart_org_path, 
                                      sep='\t', 
                                      header=0, 
                                      names=['gene_name','gene_id'])
biomart_org_gene_map = {row['gene_id']:row['gene_name'] for index, row in biomart_org_gene_map_df.iterrows()}
valid_sym_biomartorg = [biomart_org_gene_map[gid] for gid in tf.ensg if gid in biomart_org_gene_map]
sym_biomartorg_in_gtex = [s for s in valid_sym_biomartorg if s in gtex_genes]
print('#GeneSymbols from biomart.org mapped do GTEx: ' + str(len(sym_biomartorg_in_gtex)))

