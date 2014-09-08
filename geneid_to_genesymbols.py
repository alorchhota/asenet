import pandas as pd
import re
import time
#import gtf_to_genes # src: https://code.google.com/p/gtf-to-genes/

''' Show initial timestamp '''
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

'''read gencode annotation file'''

gencode_annotation_src = 'data/gencode.v19.annotation.gtf'
#gencode_annotation_src = 'data/small.annotation.gtf'
gencode_processed_annotation_dest = 'results/processed-annotations.txt'
geneid_to_genename_map_dest = 'results/geneid_gene_name.tx' 
 
'''mandatory fields in gencode data'''
ann_names = ['chr',
             'src', 
             'feature', 
             'start', 
             'end', 
             'score', 
             'strand', 
             'phase', 
             'info']

annotations = pd.read_table(gencode_annotation_src, 
                            sep='\t', 
                            comment='#', 
                            header=None, 
                            nrows=None, 
                            names=ann_names, 
                            skiprows=5)

print(len(annotations))

'''Handle extra line or carriage return at the end of a file'''
annotations = annotations.loc[0:annotations.last_valid_index(),:]


def parse_gencode_additional_info(info, attr):
    '''extract attributes for additional information of an entity'''
    regex = '([^\s]+)\s([^;]+);'
    matches = {m.group(1): m.group(2).strip('\s\t\n\r"')\
               for m in re.finditer(regex, info)}
    values = [matches[a] if a in matches else '' for a in attr] 
    return(values)

'''additional attributes in gencode data'''
attr = ['gene_id', 
        'transcript_id', 
        'gene_type', 
        'gene_status',
        'gene_name', 
        'transcript_type', 
        'transcript_status', 
        'transcript_name', 
        'exon_number', 
        'exon_id', 
        'level']

parsed_additional_info = [parse_gencode_additional_info(add_info, attr) 
                          for add_info in annotations['info']]

frame = pd.DataFrame(parsed_additional_info, columns=attr)
for a in attr:
    annotations[a] = frame[a]

'''write processed annotations to a file'''
proc_attr =  ann_names + attr
proc_attr =  [a for a in proc_attr if a != 'info']
annotations.to_csv(path_or_buf=gencode_processed_annotation_dest, columns=proc_attr, sep='\t', header=True)


''' Show end timestamp '''
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

def map_col(frame, fromCol, toCol):
    frame_by_col1 = frame.groupby(fromCol)
    frame_by_col1 = frame_by_col1.apply(lambda x: iter(x[toCol]).next())
    maps = pd.DataFrame.from_items([(fromCol, frame_by_col1.index.values), (toCol,list(frame_by_col1))])
    return(maps)

def is_col_map_unique(frame, col1, col2):
    '''checks if every item in col1 maps to exactly one item in col2'''
    col1_counts = frame.groupby(col1).apply(lambda x:len(set(x[col2])))
    return(all(col1_counts == 1))

#groupby(annotations['gene_name'], key=)
'''
if(is_col_map_unique(annotations, "gene_id", "gene_name")):
    print('OK: every gene id uniquely maps to a gene name')
else:
    print('Error: at least one gene id maps to multiple gene names')
'''

'''map gene id to gene symbol and save in a file'''
geneid_to_genename_map = map_col(annotations, 'gene_id', 'gene_name')
geneid_to_genename_map.to_csv(path_or_buf=geneid_to_genename_map_dest, sep='\t', header=True, index=False)


print('program ends!')

''' Show end timestamp '''
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
