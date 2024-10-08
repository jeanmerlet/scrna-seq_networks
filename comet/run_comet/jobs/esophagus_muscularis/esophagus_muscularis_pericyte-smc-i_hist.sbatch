#!/bin/bash

#SBATCH -A syb114
#SBATCH -J comet
#SBATCH -o jobs/esophagus_muscularis/logs/esophagus_muscularis_pericyte-smc-i_comet.%j.out
#SBATCH -e jobs/esophagus_muscularis/logs/esophagus_muscularis_pericyte-smc-i_comet.%j.err
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -p batch

executable="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/install_single_release_frontier/bin/genomics_metric"
launch_command="env OMP_PROC_BIND=spread OMP_PLACES=sockets OMP_NUM_THREADS=7 srun -N1 -n8 --cpus-per-task=7 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=closest -u"

tped_path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/esophagus_muscularis/healthy/gtex/tped/esophagus_muscularis_pericyte-smc-i/esophagus_muscularis_pericyte-smc-i_comet-mtx.bin"
comet_output_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/esophagus_muscularis/healthy/gtex/comet_out/esophagus_muscularis_pericyte-smc-i"

$launch_command $executable --num_way 2 --metric_type duo --sparse yes --num_vector 18399 \
    --num_field 14000 --all2all yes --compute_method GPU --num_proc_vector 8 --num_proc_field 1 \
    --num_proc_repl 1 --tc 1 --ccc_param 0.0 --num_tc_steps 1 --num_phase 1 --num_stage 1 --verbosity 1 \
    --checksum no --phase_min 0 --phase_max 0 --stage_min 0 --stage_max 0 --metrics_shrink 5 \
    --input_file $tped_path --histograms_file $comet_output_dir/hist.tsv --threshold 4

echo "finished esophagus_muscularis_pericyte-smc-i in esophagus_muscularis"