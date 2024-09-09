import subprocess
import os, re
import pandas as pd
import numpy as np


data_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human'
jobs_dir = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/run_comet/jobs'


# settings
histogram = False
threshold = True
pct = 0.20


# write a jobscript to file
# this is formatted for frontier
def submit_job(tped_path, jobs_dir, out_dir, celltype, tissue, num_vectors, num_fields, histogram, threshold, ccc=False):
    if ccc:
        ccc_txt = ''
    else:
        ccc_txt = '--ccc_param 0.0 '
    jobscript = (
        '#!/bin/bash\n'
        '\n'
        '#SBATCH -A syb114\n'
        '#SBATCH -J comet\n'
        f'#SBATCH -o jobs/{tissue}/logs/{celltype}_comet.%j.out\n'
        f'#SBATCH -e jobs/{tissue}/logs/{celltype}_comet.%j.err\n'
        '#SBATCH -t 00:30:00\n'
        '#SBATCH -N 1\n'
        '#SBATCH -p batch\n'
        '\n'
        'executable="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/install_single_release_frontier/bin/genomics_metric"\n'
        'launch_command="env OMP_PROC_BIND=spread OMP_PLACES=sockets OMP_NUM_THREADS=7'
        ' srun -N1 -n8 --cpus-per-task=7 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=closest -u"\n'
        '\n'
        f'tped_path="{tped_path}"\n'
        f'comet_output_dir="{out_dir}"\n'
        '\n'
        f'$launch_command $executable --num_way 2 --metric_type duo --sparse yes --num_vector {num_vectors} \\\n'
        f'    --num_field {num_fields} --all2all yes --compute_method GPU --num_proc_vector 8 --num_proc_field 1 \\\n'
        f'    --num_proc_repl 1 --tc 1 {ccc_txt}--num_tc_steps 1 --num_phase 1 --num_stage 1 --verbosity 1 \\\n'
        '    --checksum no --phase_min 0 --phase_max 0 --stage_min 0 --stage_max 0 --metrics_shrink 5 \\\n'
        '    --input_file $tped_path'
    )
    jobs_dir = os.path.join(jobs_dir, tissue)
    logs_dir = os.path.join(jobs_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    if histogram:
        out_str = '--histograms_file $comet_output_dir/hist.tsv'
        out_str += ' --threshold 4'
        job_path = os.path.join(jobs_dir, celltype + '_hist.sbatch')
    else:
        if ccc:
            out_str = '--output_file_stub $comet_output_dir/out_ccc'
        else:
            out_str = '--output_file_stub $comet_output_dir/out'
        if threshold:
            hist_path = os.path.join(out_dir, 'hist.tsv')
            cutoff = get_cutoff(hist_path, pct)
            out_str += f' --threshold 4,4,4,{cutoff}'
            job_path = os.path.join(jobs_dir, celltype + '_cutoff-' + str(cutoff) + '.sbatch')
        else:
            job_path = os.path.join(jobs_dir, celltype + '.sbatch')
    jobscript = ' '.join((jobscript, out_str))
    jobscript = '\n\n'.join((jobscript, f'echo "finished {celltype} in {tissue}"'))
    with open(job_path, 'w') as job_file:
        job_file.write(jobscript)
    subprocess.run(['sbatch', job_path])


# count number of input vectors and their length
def count_dims(tped_path):
    tped_path = tped_path[:-3] + 'tped'
    num_vectors = subprocess.run(['wc','-l', tped_path], capture_output=True, text=True).stdout.split(' ')[0]
    p1 = subprocess.Popen(['head','-n1', tped_path], stdout=subprocess.PIPE)
    p2 = subprocess.run(['wc', '-w'], stdin=p1.stdout,capture_output=True, text=True)
    p1.wait()
    num_fields = str(int(p2.stdout.split(' ')[0])-4)
    return num_vectors, num_fields


# determine threshold based on histogram
# pct is given as a proportion
def get_cutoff(hist_path, pct):
    hist = pd.read_csv(hist_path, sep='\t')
    ll_hh = hist.iloc[:, 4].values
    total_edges = np.sum(ll_hh)
    target_num_edges = int(pct * total_edges)
    for idx, num_edges in enumerate(ll_hh):
        total_edges -= num_edges
        if total_edges < target_num_edges:
            idx -= 1
            break
    cutoff = hist.iloc[idx, 0]
    return cutoff


def comet_already_run(out_dir):
    for f in os.listdir(out_dir):
        if '.bin' in f:
            return True
    return False


# get a list of absolute tped paths
input_paths = []
for r, d, f in os.walk(data_dir):
    for tped in f:
        if '-mtx.bin' in tped:
            input_paths.append(os.path.join(r, tped))

input_paths.sort()

# write and submit jobs by celltype
start_i = 200
end_i = 300
for i, path in enumerate(input_paths):
    if i < start_i: continue
    if i == end_i: break
    prefix = re.search('^(.*)/tped/.*$', path).groups()[0]
    celltype = re.search('tped/(.*)/', path).groups()[0]
    tissue = re.search('human/(.*)/healthy', path).groups()[0]
    root_out_dir = os.path.join(prefix, 'comet_out')
    out_dir = os.path.join(root_out_dir, celltype)
    os.makedirs(root_out_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    num_vectors, num_fields = count_dims(path)
    if not comet_already_run(out_dir):
        submit_job(path, jobs_dir, out_dir, celltype, tissue, num_vectors, num_fields, histogram, threshold, ccc=True)
        print(f'submitted {celltype} in {tissue}', flush=True)
