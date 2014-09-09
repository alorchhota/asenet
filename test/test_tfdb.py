import tfdb
import os

os.chdir('/home/ashis/work/github/asenet/')


src = 'data/tf-nrg2538-s3.txt'
tf = tfdb.TFDB(src)
print(len(tf.sym))
valid_sym = [s!='NA' for s in tf.sym]
print(sum(valid_sym))
