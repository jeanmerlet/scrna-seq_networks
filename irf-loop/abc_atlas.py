
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import anndata

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
from abc_atlas_access.abc_atlas_cache.anndata_utils import get_gene_data

download_base = Path('/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/brain/healthy/allen_brain/abc_atlas')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

# list current release
abc_cache.current_manifest
# 'releases/20240831/manifest.json'

# list previous releases
abc_cache.list_manifest_file_names

#abc_cache.load_manifest('releases/20230630/manifest.json')
#print("old manifest loaded:", abc_cache.current_manifest)

# Return to the latest manifest
#abc_cache.load_latest_manifest()
#print("after latest manifest loaded:", abc_cache.current_manifest)

# list available directories
abc_cache.list_directories

# list data files available in a directory
abc_cache.list_data_files('WHB-10Xv3')
# ['WHB-10Xv3-Neurons/log2', 'WHB-10Xv3-Neurons/raw', 'WHB-10Xv3-Nonneurons/log2', 'WHB-10Xv3-Nonneurons/raw']

# we're going to want to download both raw neuron and non-neuron databases
# liat metadata files available in a directory
abc_cache.list_metadata_files('WHB-10Xv3')
# ['anatomical_division_structure_map', 'cell_metadata', 'donor', 'example_genes_all_cells_expression', 'gene', 'region_of_interest_structure_map']

# get size of directory
abc_cache.get_directory_data_size('WHB-10Xv3')
# '70.03 GB'

##### expression data

# download gene expression matrix
abc_cache.get_data_path(directory='WHB-10Xv3',file_name='WHB-10Xv3-Neurons/raw')
abc_cache.get_data_path(directory='WHB-10Xv3',file_name='WHB-10Xv3-Nonneurons/raw')

##### metadata

# download cell metadata
cell = abc_cache.get_metadata_dataframe(directory='WHB-10Xv3', file_name='cell_metadata').set_index('cell_label')
cell.head
cell.columns

# download donor metadata
donor = abc_cache.get_metadata_dataframe(directory='WHB-10Xv3', file_name='donor')
donor.head
donor.columns

# download anatomical metadata
anatom = abc_cache.get_metadata_dataframe(directory='WHB-10Xv3', file_name='anatomical_division_structure_map')
anatom.head
anatom.columns

# download region metadata
region = abc_cache.get_metadata_dataframe(directory='WHB-10Xv3', file_name='region_of_interest_structure_map')
region.head
region.columns

# download gene metadata
gene = abc_cache.get_metadata_dataframe(directory='WHB-10Xv3', file_name='gene').set_index('gene_identifier')
gene.head
gene.columns

# download cluster annotations
abc_cache.list_metadata_files('WHB-taxonomy')
cluster = abc_cache.get_metadata_dataframe(directory='WHB-taxonomy', file_name='cluster')
cluster_anno = abc_cache.get_metadata_dataframe(directory='WHB-taxonomy', file_name='cluster_annotation_term')
cluster_anno_set = abc_cache.get_metadata_dataframe(directory='WHB-taxonomy', file_name='cluster_annotation_term_set')
cluster_mem = abc_cache.get_metadata_dataframe(directory='WHB-taxonomy', file_name='cluster_to_cluster_annotation_membership')




