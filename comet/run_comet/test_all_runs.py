import os, re


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'


# get a list of absolute tped paths
input_paths = []
for r, d, f in os.walk(data_dir):
    for tped in f:
        if '-mtx.bin' in tped:
            input_paths.append(os.path.join(r, tped))

input_paths.sort()
print(len(input_paths))
raise SystemExit()


def comet_already_run(out_dir):
    for f in os.listdir(out_dir):
        if '.bin' in f:
            return True
    return False



# write and submit jobs by celltype
for i, path in enumerate(input_paths):
    prefix = re.search('^(.*)/tped/.*$', path).groups()[0]
    celltype = re.search('tped/(.*)/', path).groups()[0]
    #if celltype == 'schwann-cell-i': continue
    tissue = re.search('human/(.*)/healthy', path).groups()[0]
    root_out_dir = os.path.join(prefix, 'comet_out')
    out_dir = os.path.join(root_out_dir, celltype)
    if not comet_already_run(out_dir):
        print(f'{celltype} in {tissue} not run', flush=True)
