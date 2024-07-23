import loompy
import os

hypo_looms_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/nemo/raw/linnarson/hypo_looms'
hypo_loom_paths = [os.path.join(hypo_looms_dir, f) for f in os.listdir(hypo_looms_dir)]
hypo_out_path = os.path.join(hypo_looms_dir, 'all_hypothalamus.loom')

loompy.combine(files=hypo_loom_paths, output_file=hypo_out_path, key='Accession', batch_size=500)
#loompy.combine(files=hypo_loom_paths, output_file=hypo_out_path, batch_size=10000)
