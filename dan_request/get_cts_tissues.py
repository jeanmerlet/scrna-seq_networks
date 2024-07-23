import os, re
import pandas as pd

tissue_dirs = ['/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_mtx',
               '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/liver/healthy/vcir/comet_mtx',
               '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/lung/healthy/gtex/comet_mtx']
tissue_names = ['brain', 'liver', 'lung']

cts_per_tissue = {}
for i, tissue_dir in enumerate(tissue_dirs):
    ct_filenames = os.listdir(tissue_dir)
    tissue_name = tissue_names[i] 
    for ct_filename in ct_filenames:
        if 'brain' in tissue_dir:
            if 'ba8' in ct_filename or 'ba9' in ct_filename or 'ba10' in ct_filename: continue
            ct = re.search('pfc_(.*)_comet-mtx', ct_filename).groups()[0]
        else:
            ct = ct_filename.split('_')[1]
        if tissue_name in cts_per_tissue.keys():
            cts_per_tissue[tissue_name].append(ct)
        else:
            cts_per_tissue[tissue_name] = [ct]

tissue_col = []
ct_col = []

for tissue, cts in cts_per_tissue.items():
    for ct in cts:
        if tissue == 'brain':
            tissue_col.append('brain-pfc')
        else:
            tissue_col.append(tissue)
        ct_col.append(ct)

df = pd.DataFrame(data=[tissue_col, ct_col]).T


out_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/meta/brain-liver-lung_cts.tsv'
df.to_csv(out_path, sep='\t', header=None, index=None)

 


