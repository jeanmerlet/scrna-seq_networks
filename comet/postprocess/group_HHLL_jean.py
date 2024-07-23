import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir')
parser.add_argument('-r', '--run', action='store_true')
args = parser.parse_args()


data_paths = []
for r, d, f in os.walk(args.data_dir):
    for txt in f:
        if 'comet_out' in r and '.txt' in txt:
            if 'archive' not in r:
                data_paths.append(os.path.join(r, txt))
data_paths.sort()


def store_edgeweight(name, HH_vals, LL_vals, duo0, duo1, edgeweight):
    if duo0 == '1' and duo1 == '1':
        HH_vals[name] = edgeweight
    elif duo0 == '0' and duo1 == '0':
        LL_vals[name] = edgeweight

def combine_edgeweights(name, HH_vals, LL_vals, HHLL_vals):
    if name in HH_vals and name in LL_vals:
        combined_edgeweight = HH_vals[name] + LL_vals[name]
        HHLL_vals[name] = combined_edgeweight

def convert_to_edge_list_tsv(in_path):
    HH_vals, LL_vals, HHLL_vals = {}, {}, {}
    out_dir, file_name = os.path.split(path)
    out_path = os.path.join(out_dir, 'HHLL_' + file_name[:-4] + '.tsv')
    with open(in_path, 'rt') as in_file:
        for line in in_file:
            _, duo0, _, duo1, id0, id1, edgeweight = line.strip().split(' ')
            edgeweight = float(edgeweight)
            if duo0 == duo1:
                name = [id0[:-2], id1[:-2]]
                name.sort()
                name = tuple(name)
                store_edgeweight(name, HH_vals, LL_vals, duo0, duo1, edgeweight)
                combine_edgeweights(name, HH_vals, LL_vals, HHLL_vals)
    with open(out_path, 'wt') as out_file:
        for name, edgeweight in HHLL_vals.items():
            out_file.write('\t'.join([name[0], name[1], str(round(edgeweight, 6))]) + '\n')
                    

if not args.run:
    print(f'txt files to convert to combined tsv: {len(data_paths)}')
else:                    
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    for i, path in enumerate(data_paths):
        # match to ranks
        if i != rank: continue
        convert_to_edge_list_tsv(path)
        print(f'finished {path}')
