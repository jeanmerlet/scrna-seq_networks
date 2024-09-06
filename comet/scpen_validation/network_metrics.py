# import packages
import networkx as nx
import pandas as pd
import numpy as np

# function to create a network from pandas data frame
def create_network(file):
    data = pd.read_csv(file,sep = "\t",header = None)
    G = nx.Graph()
    edges = list(data.itertuples(index = False,name = None))
    G.add_weighted_edges_from(edges)
    return(G)

# function to subset a graph
def node_sublist(network,indices):
    nodes = list(network.nodes)
    nodes = [nodes[index] for index in indices]
    return(nodes)

# function to calculate connectivity
# TODO(?): sum of edge weights later
def get_connectivity(network,indices):
    nodes = node_sublist(network,indices)
    connectivity_result = {}
    for node in nodes:
         node_subgraph = network.edge_subgraph(network.edges(node))
         connectivity = nx.node_connectivity(node_subgraph)
         connectivity_result[node] = connectivity
    return(connectivity_result)
    
# function to calculate clustering coefficient
def get_clustering_coeff(network,indices):
    nodes = node_sublist(network,indices)
    clustering_coeff_result = {}
    for node in nodes:
         node_subgraph = network.edge_subgraph(network.edges(node))
         nodes_subgraph = node_subgraph.nodes
         subgraph = network.subgraph(nodes_subgraph)
         clustering_coeff_result.update(nx.clustering(subgraph))
    return(clustering_coeff_result)

def get_metrics(file,indices,out_path):
    network = create_network(file)
    connectivity_result = get_connectivity(network,indices)
    clustering_coeff_result = get_clustering_coeff(network,indices)
    index = []
    data = np.zeros((len(connectivity_result.keys()),2),dtype = int)
    for i,key in enumerate(connectivity_result.keys()):
        index.append(key)
        data[i,0] = connectivity_result[key]
        data[i,1] = clustering_coeff_result[key]
    metrics = pd.DataFrame(data = data,index = index,columns = ['connectivity','clustering_coefficient'])
    metrics.to_csv(out_path, sep='\t')
