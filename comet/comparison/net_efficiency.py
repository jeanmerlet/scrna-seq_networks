import numpy as np
import pandas as pd
import os, re
import networkx as nx


irf_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/irf-loop_out/immune-dc-macrophage/heart_immune-dc-macrophage_irf-loop-mtx_no_correlated_data_top_0.01_pct_correlates_readded.tsv'
comet_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_out/myocyte-cardiac/top_1_combined.tsv'

irf_df = pd.read_csv(irf_path, sep='\t', header=None)
comet_df = pd.read_csv(comet_path, sep='\t', header=None)
gene_count=[]
def make_network(in_file):
    gene_a = list(in_file.iloc[:,0])
    gene_b = list(in_file.iloc[:,1])
    weight = list(in_file.iloc[:,2])
    G = nx.Graph()
    for gene in gene_a:
        if gene not in gene_count:
            gene_count.append(gene)
    for gene in gene_b:
        if gene not in gene_count:
            gene_count.append(gene)

    print(len(gene_count))

    for i, weight in enumerate(weight):
        G.add_edge(gene_a[i], gene_b[i], weight = weight)
    return G

comet_net = make_network(comet_df)
#irf_net = make_network(irf_df)


def get_efficiency(net):
    efficiency_list= []
    for i in nx.nodes(net):
        for j in list(nx.nodes(net)):
            print(j)
            if i == j or j in list(nx.nodes(net))[0:i]:
                continue
            eff = 1 / nx.shortest_path_length(net, i, j,weight='weight')
            efficiency_list.append(eff)

    return efficiency_list
#irf_efficiency = get_efficiency(irf_net)
#print(irf_efficiency)
