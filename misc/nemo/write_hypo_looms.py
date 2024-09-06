import loompy
import numpy as np
import scipy.sparse as sparse
import os

path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/nemo/raw/linnarson/adult_human_20221007.loom'

with loompy.connect(path) as ds:
    coarse = ds.ca['ROIGroupCoarse']
    hypo_idx = coarse == 'Hypothalamus'
    all_bcs = ds.ca['CellID']
    hypo_bcs = all_bcs[hypo_idx]

hypo_sample_ids = [bc[:8] for bc in hypo_bcs]
uniq_hypo_sample_ids, hypo_sample_counts = np.unique(hypo_sample_ids, return_counts=True)

hypo_bcs_by_sample = {}
for hypo_sample_id in uniq_hypo_sample_ids:
    hypo_sample_bcs = np.array([True if bc[:8] == hypo_sample_id else False for bc in all_bcs])
    hypo_bcs_by_sample[hypo_sample_id] = hypo_sample_bcs

with loompy.connect(path) as ds:
    ra = {}
    for key, value in ds.ra.items():
        ra[key] = value

hypo_loom_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/nemo/raw/linnarson/hypo_looms'
i = 0
with loompy.connect(path) as ds:
    for sample, idx in hypo_bcs_by_sample.items():
        i += 1
        if i < 7: continue
        print(sample)
        out_path = os.path.join(hypo_loom_dir, str(sample) + '.loom')
        ca = {}
        for key, value in ds.ca.items():
            ca[key] = value[idx]
        mtx = ds[:, idx]
        loompy.create(out_path, mtx, ra, ca)
