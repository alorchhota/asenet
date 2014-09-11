import gtf
import pandas as pd
import os

home_dir = '/home/asaha6/github/asenet'

# set working directory
os.chdir(home_dir)

gencode_annotation_src = 'data/gencode.v19.annotation.gtf'
gr = gtf.GencodeGTFReader(gencode_annotation_src, False)
print(len(gr.data))
#print(gr.data.ix[0:4])
processed_src = 'results/gencode.v19.annotation.processed.gtf'
gr.save(processed_src)


#processed_src = '/home/ashis/work/github/asenet/data/gencode.v19.annotation.processed.gtf' 
#gr2 = gtf.GencodeGTFReader(processed_src, processed=True)

#dic = gr2.map_geneid_to_genename()
#print(dic['ENSG00000223972'])

#df = pd.DataFrame.from_items([('gene_id',dic.keys()),('gene_name',dic.values())])
#df.to_csv(path_or_buf='results/genemap.txt', sep='\t', header=True)

