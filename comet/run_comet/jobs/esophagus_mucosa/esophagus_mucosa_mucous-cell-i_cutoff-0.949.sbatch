#!/bin/bash

#SBATCH -A syb111
#SBATCH -J comet
#SBATCH -o jobs/esophagus_mucosa/logs/esophagus_mucosa_mucous-cell-i_comet.%j.out
#SBATCH -e jobs/esophagus_mucosa/logs/esophagus_mucosa_mucous-cell-i_comet.%j.err
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -p batch

executable="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/install_single_release_frontier/bin/genomics_metric"
launch_command="env OMP_PROC_BIND=spread OMP_PLACES=sockets OMP_NUM_THREADS=7 srun -N1 -n8 --cpus-per-task=7 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=closest -u"

tped_path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/esophagus_mucosa/healthy/gtex/tped/esophagus_mucosa_mucous-cell-i/0.000000001_var_esophagus_mucosa_mucous-cell-i_comet-mtx.bin"
comet_output_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/esophagus_mucosa/healthy/gtex/comet_out/esophagus_mucosa_mucous-cell-i"

$launch_command $executable --num_way 2 --metric_type duo --sparse yes --num_vector 23461 \
    --num_field 6160 --all2all yes --compute_method GPU --num_proc_vector 8 --num_proc_field 1 \
    --num_proc_repl 1 --tc 1 --num_tc_steps 1 --num_phase 1 --num_stage 1 --verbosity 1 \
    --checksum no --phase_min 0 --phase_max 0 --stage_min 0 --stage_max 0 --metrics_shrink 5 \
    --input_file $tped_path --output_file_stub $comet_output_dir/out_ccc --threshold 4,4,4,0.949

echo "finished esophagus_mucosa_mucous-cell-i in esophagus_mucosa"