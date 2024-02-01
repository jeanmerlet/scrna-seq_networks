module load cray-python/3.10.10
#normalize scores to 0-1 divide by largest num

data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/heart/healthy/gtex/comet_out/adipocyte'

data_paths = [os.path.join(data_dir, path) for path in os.listdir(data_dir) if 'top' in path]
data_paths.sort()


def normalize_scores(in_path):
    _, file_name = os.path.split(path)
    out_path = os.path.join(out_dir, 'normalized_' + file_name)
    with open(out_path, 'wt') as out_file:
        with open(in_path, 'rt') as in_file:
            for line in in_file:
                id0, id1, edge_weight = line.strip('\n').split(' ')
                edge_weight = float(edge_weight)
                name = [id0, id1]
                name.sort()
                name = tuple(name)
                norm_weights = {}
                norm_weights[name] = edge_weight

            for i, key in enumerate(norm_weights):
                out_file.write('\t'.join((key[0], key[1], str(round(HHLL_vals[key],6)))) + '\n')



comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

for i, path in enumerate(data_paths):
    # distribute across ranks
    if i % size != rank: continue
    normalize_scores(path)
