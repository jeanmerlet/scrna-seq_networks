import re, os


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'


file_list = []
for r, d, f in os.walk(data_dir):
    for mtx in f:
        if 'irf-loop-mtx.tsv' in mtx:
            if 'archive' not in r:
                if 'brain' not in r:
                    if '.err' not in mtx:
                        if '.out' not in mtx:
                            if 'Archive' not in r:
                                file_list.append(os.path.join(r, mtx))
file_list.sort()
print(len(file_list))


out_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/for_outlier_removal/matrices.txt'

with open(out_path, 'w') as out_file:
    for f in file_list:
        out_file.write(f)
        out_file.write('\n')
