import numpy as np
import pandas as pd
import os, re

comet_mtx_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_mtx/heart_immune-dc-macrophage_comet-mtx.tsv'

irf_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/irf-loop_out/immune-dc-macrophage/heart_immune-dc-macrophage_irf-loop-mtx_no_correlated_data_top_0.01_pct_correlates_readded.tsv'
#comet_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_out/immune-dc-macrophage/top_1_combined_out.tsv'
comet_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_out/immune-dc-macrophage/edge_list/top-1_HH+LL_combined_out.tsv'


irf_df = pd.read_csv(irf_path, sep='\t', header=None)
comet_df = pd.read_csv(comet_path, sep='\t', header=None)


irf_gene_pairs = []
comet_gene_pairs = []

def read_gene_pairs(irf_path, comet_path, bidirectional=False):
    with open(irf_path, 'rt') as irf_file:
        for irf_line in irf_file:
            irf_gene_a, irf_gene_b, irf_weight = irf_line.strip().split('\t')
            gene_pair = ' '.join((irf_gene_a, irf_gene_b))
            irf_gene_pairs.append(gene_pair)
    with open(comet_path, 'rt') as comet_file:
        for comet_line in comet_file:
            comet_gene_a, comet_gene_b, comet_weight = comet_line.strip().split('\t')
            gene_pair = re.sub('_T', '', re.sub('_A', '', ' '.join((comet_gene_a, comet_gene_b))))
            comet_gene_pairs.append(gene_pair)
            if bidirectional:
                rev_gene_pair = re.sub('_T', '', re.sub('_A', '', ' '.join((comet_gene_b, comet_gene_a))))
            comet_gene_pairs.append(gene_pair)


read_gene_pairs(irf_path, comet_path, bidirectional=True)
print('gene pairs read')

irf_gene_pairs.sort()
comet_gene_pairs.sort()

irf_gene_pairs = np.array(irf_gene_pairs)
comet_gene_pairs = np.array(comet_gene_pairs)

print(irf_gene_pairs[:20])
print(comet_gene_pairs[:20])
print(len(irf_gene_pairs))
print(len(comet_gene_pairs))

matching_gene_paris = np.intersect1d(irf_gene_pairs, comet_gene_pairs)
num_matching_gene_pairs = len(np.intersect1d(comet_gene_pairs, irf_gene_pairs))

print(num_matching_gene_pairs)


irf_genes = []
for gene_pair in irf_gene_pairs:
    gene_a, gene_b = gene_pair.split(' ')
    irf_genes.append(gene_a)
    irf_genes.append(gene_b)

comet_genes = []
for gene_pair in comet_gene_pairs:
    gene_a, gene_b = gene_pair.split(' ')
    comet_genes.append(gene_a)
    comet_genes.append(gene_b)

print(f'# uniq irf genes: {len(np.unique(irf_genes))}')
print(f'# uniq comet genes: {len(np.unique(comet_genes))}')
print(len(np.intersect1d(comet_genes, irf_genes)))

print('reading comet mtx')
comet_mtx = pd.read_csv(comet_mtx_path, sep='\t')
comet_mtx.index = [x.strip('"') for x in list(comet_mtx.index.values)]
print('read comet mtx')

for gene in np.unique(comet_genes)[:20]:
    print(f'{gene} total: {np.sum(comet_mtx.loc[gene, :].values)}')
