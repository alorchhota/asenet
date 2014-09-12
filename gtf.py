import pandas as pd
import re
import numpy as np

class GencodeGTFReader:
    
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
    
    proc_attr =  [a for a in (ann_names + attr) if a != 'info']


    def parse_gencode_additional_info(self, info, attr):
            '''extract attributes for additional information of an entity'''
            regex = '([^\s]+)\s([^;]+);'
            matches = {m.group(1): m.group(2).strip('\s\t\n\r"')\
                       for m in re.finditer(regex, info)}
            values = [matches[a] if a in matches else '' for a in attr] 
            return(values)
    
    def __init__(self, src, processed=False):
        '''
        initialize gtf data from source file.
        src:         source file url
        processed:   if processed=False, loads from raw file
                     otherwise load from processed file.
        '''
        
        self.data = None    # contains processed data
        if processed is False:
            self.populate_from_raw_file(src)
        else:
            self.populate_from_processed_file(src)
        
    def populate_from_raw_file(self, src):
        annotations = pd.read_table(src, 
                            sep='\t', 
                            comment='#', 
                            header=None, 
                            nrows=None, 
                            names= GencodeGTFReader.ann_names,
                            skiprows=5)
         
        '''Handle extra line or carriage return at the end of a file'''
        if len(annotations) > annotations.last_valid_index()+1:
            annotations = annotations.loc[0:annotations.last_valid_index(),:]
        
        
 
        parsed_additional_info = [self.parse_gencode_additional_info(add_info, GencodeGTFReader.attr)
                                  for add_info in annotations['info']]

        frame = pd.DataFrame(parsed_additional_info, columns=GencodeGTFReader.attr)
        for a in GencodeGTFReader.attr:
            annotations[a] = frame[a]
            
        self.data = annotations[GencodeGTFReader.proc_attr]
        

    def populate_from_processed_file(self, src):
        annotations = pd.read_table(src, 
                            sep='\t', 
                            comment='#', 
                            header=0, 
                            nrows=None)
         
        '''Handle extra line or carriage return at the end of a file'''
        if len(annotations) > annotations.last_valid_index()+1:
            annotations = annotations.loc[0:annotations.last_valid_index(),:]
        
        self.data = annotations
        
    def save(self, dest, ext='tsv'):
        '''save the annotation data in a file'''
        sep = '\t'
        if ext=='csv':
            sep = ','
        
        self.data.to_csv(path_or_buf=dest, columns=GencodeGTFReader.proc_attr, sep=sep, header=True)
        
    def map_col(self, fromCol, toCol):
        ''' map from one column in data to another (only one value) '''
        frame_by_col1 = self.data.groupby(fromCol)
        frame_by_col1 = frame_by_col1.apply(lambda x: iter(x[toCol]).next())
        
        col1_values = frame_by_col1.index.values
        col2_values = list(frame_by_col1)
        
        maps = {col1_values[i]:col2_values[i] for i in range(len(col1_values))}
        return(maps)

    def map_geneid_to_genename(self):
        fromCol, toCol = 'gene_id', 'gene_name'
        frame_by_col1 = self.data.groupby(fromCol)
        frame_by_col1 = frame_by_col1.apply(lambda x: iter(x[toCol]).next())
        
        col1_values = frame_by_col1.index.values
        col2_values = list(frame_by_col1)
        
        maps = {col1_values[i].split('.')[0]:col2_values[i] for i in range(len(col1_values))}
        return(maps)
    
    def get_gene_name_from_locus(self, n_chr, pos):
        '''
        get gene name of a locus.
        If multiple genes found, take an exonic gene.
        '''
        
        matched_chr = self.data.loc[:,'chr'] == ('chr'+str(n_chr))
        matched_start = self.data.loc[:,'start'] <= pos
        matched_end = self.data.loc[:,'end'] >= pos
        
        matched = matched_chr & matched_start & matched_end
        matched_index = np.where(matched)
        matched_entities = self.data.loc[self.data.index.values[matched_index], ['gene_name', 'feature']]
        
        gene_names = set(matched_entities.loc[:,'gene_name'])
        
        if len(gene_names) == 0:
            return '-'
        elif len(gene_names) == 1:
            return gene_names.pop()
        else:
            is_exon = matched_entities.loc[:,'feature'] == 'exon'
            if sum(is_exon) > 0:
                exon_index = np.where(is_exon)
                exon_gene = matched_entities.iloc[exon_index[0][0], 0]
                return exon_gene
            return gene_names.pop()
        
    