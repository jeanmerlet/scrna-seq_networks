import pandas as pd
import os

dir_path = "/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_out/pfc_glial-non-neuron_micro/networks/rwr"
out_auroc_path = "/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_out/pfc_glial-non-neuron_micro/networks/summary_auroc_doubles.tsv"

def remove_network_name(network_file): 
    prefix_removed = network_file.split("top_")[1]
    suffix_removed = prefix_removed.split("_ensembl")[0]
    return suffix_removed


network_paths = [os.path.join(dir_path, path) for path in os.listdir(dir_path)]
total_df = {}
for network_path in network_paths:
    cv_outputs = os.listdir(network_path)
    output_files = list(filter(lambda x: "summary" in x, cv_outputs))    
    series_data = {}
    for summary_file in output_files: 
        summary_df = pd.read_csv(f'{network_path}/{summary_file}', sep='\t')
        auroc_data = summary_df[summary_df['measure'] == 'AUROC']
        mean_auroc = auroc_data[auroc_data['fold'] == 'meanrank']['value'].mean()
        geneset = auroc_data['geneset'].values[0]
        series_data[geneset] = mean_auroc
    #network_name = remove_network_name(network_type)
    _, network_name = os.path.split(network_path)
    total_df[network_name] = series_data

total_df = pd.DataFrame(total_df)
total_df.to_csv(out_auroc_path, sep='\t') 
