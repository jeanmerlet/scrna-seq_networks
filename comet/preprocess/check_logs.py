import os, re


log_path = '/lustre/orion/syb111/proj-shared/Projects/scrna-seq/code/comet/preprocess/logs/comet_mtx_to_tped.2370421.out'

finished = []
with open(log_path) as logfile:
    for line in logfile:
        if 'Finished' in line:
            #print(line.strip())
            name = re.search('from (.*) \(', line.strip()).groups()[0]
            finished.append(name)


unfinished = []
with open(log_path) as logfile:
    for line in logfile:
        if 'Finished' not in line:
            #print(line.strip())
            name = re.search('from (.*) \(', line.strip()).groups()[0]
            if name not in finished:
                unfinished.append(name)

for name in unfinished:
    print(name)
