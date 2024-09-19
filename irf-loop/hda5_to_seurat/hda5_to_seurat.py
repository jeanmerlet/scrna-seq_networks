
#example script
#to run:
#mkdir <out_dir>
#python name_of_this_script.py <adata_path> <out_dir>
#gzip <out_dir>/*

import scanpy as sc
from scipy import io
import sys

adata = sc.read_h5ad(sys.argv[1])
out_dir = sys.argv[2]

#adata = adata.raw.to_adata() #only if adata has RAW saved and thats what you want!!

with open(out_dir + '/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')

with open(out_dir + '/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')

io.mmwrite(out_dir +'/matrix.mtx', adata.X.T)

adata.obs.to_csv(sys.argv[1] + '.metadata.csv')

