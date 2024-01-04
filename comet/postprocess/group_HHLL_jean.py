from mpi4py import MPI
import os, re


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'

for r, d, f in os.walk(data_dir):
    for comet_outfile in f:
        if 

#out1, out2 = 'T', 'T'
out_type = 'HH+LL'

HHLL_vals, HHLL_counts,original_weights = {},{},{}

data_paths = [os.path.join(data_dir, path) for path in os.listdir(data_dir) if '.txt' in path]
data_paths.sort()




def check_correct_duo_type(duo0,duo1):
    return((duo0 == 0 and duo1 == 0) or (duo0 == 1 and duo1 == 1))

def convert_to_edge_list_tsv(in_path):
    _, file_name = os.path.split(path)
    out_path = os.path.join(out_dir, out_type + '_' + file_name[:-4] + '.tsv')
    with open(out_path, 'wt') as out_file:
        with open(in_path, 'rt') as in_file:
            for line in in_file:
                _, duo0, _, duo1, id0, id1, edge_weight = line.strip('\n').split(' ')
                duo0, duo1, edge_weight = int(duo0), int(duo1), float(edge_weight)
                if check_correct_duo_type(duo0,duo1):

                    name = [id0[:-2], id1[:-2]]
                    name.sort()
                    name = tuple(name)
                    if duo0 == 0 and duo1==0:
                        LL_vals[name] = edge_weight
                    if duo0 == 1 and duo1 ==1:
                        HH_vales[name] = edge_weight

                    if name not in HHLL_vals:
                        HHLL_vals[name] = edge_weight
                        HHLL_counts[name] = 1
                    elif HHLL_counts[name] == 1:
                        HHLL_vals[name] += edge_weight 
                        HHLL_counts[name] += 1
                    else:
                        print(name)
                        break
            
            for i, key in enumerate(HHLL_vals):
                print(key[0]+' '+key[1]+ '\n'+'HH: '+ HH_vals[key] + ' LL: '+LL_vals[key]+' HHLL: '+HHLL_vals[key]
                out_file.write('\t'.join((key[0], key[1], str(round(HHLL_vals[key],6)))) + '\n')
                    
                    
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, path in enumerate(data_paths):
    # distribute across ranks
    if i % size != rank: continue
    convert_to_edge_list_tsv(path)
