import pandas as pd

class TFDB:
    na_value = 'NA'     # this string will be used for na values
    def readNrg2538Db(self, src):
        tf_table = pd.read_table(src, sep='\t', header=0, skiprows=11, na_values='Unknown')
        if(len(tf_table) > tf_table.last_valid_index()+1):
            tf_table = tf_table.ix[0:tf_table.last_valid_index()]
        tf_table = tf_table.fillna(TFDB.na_value)
        valid_indexes = [c=='a' or c=='b' or c=='other' for c in tf_table['Class']]
        self.ensg = tf_table['Ensembl ID'][valid_indexes]
        self.sym = tf_table['HGNC symbol'][valid_indexes]
        
    def __init__(self, src, db_type=0):
        '''initialize transcriptor factor database from source file.
        src:    database file path 
        type:   0 = nature nrg2538 database. 
                Other values will be defined later if required.'''
        
        self.sym = []   # contains TF hgnc symbols 
        self.ensg = []  # contains TF ensembl gene ids
        
        if(db_type==0):
            self.readNrg2538Db(src)
        
    
        