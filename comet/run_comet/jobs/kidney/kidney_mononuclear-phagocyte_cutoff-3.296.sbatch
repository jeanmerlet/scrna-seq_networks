#!/bin/bash

#SBATCH -A syb114
#SBATCH -J comet
#SBATCH -o jobs/kidney/logs/kidney_mononuclear-phagocyte_comet.%j.out
#SBATCH -e jobs/kidney/logs/kidney_mononuclear-phagocyte_comet.%j.err
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -p batch

executable="/lustre/orion/syb111/proj-shared/Personal/jmerlet/tools/comet/install_single_release_frontier/bin/genomics_metric"
launch_command="env OMP_PROC_BIND=spread OMP_PLACES=sockets OMP_NUM_THREADS=7 srun -N1 -n8 --cpus-per-task=7 --ntasks-per-node=8 --gpus-per-task=1 --gpu-bind=closest -u"

tped_path="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/kidney/healthy/kpmp/tped/kidney_mononuclear-phagocyte/kidney_mononuclear-phagocyte_comet-mtx.bin"
comet_output_dir="/lustre/orion/syb111/proj-shared/Projects/scrna-seq/data/human/kidney/healthy/kpmp/comet_out/kidney_mononuclear-phagocyte"

$launch_command $executable --num_way 2 --metric_type duo --sparse yes --num_vector 22865 \
    --num_field 21840 --all2all yes --compute_method GPU --num_proc_vector 8 --num_proc_field 1 \
    --num_proc_repl 1 --tc 1 --ccc_param 0.0 --num_tc_steps 1 --num_phase 1 --num_stage 1 --verbosity 1 \
    --checksum no --phase_min 0 --phase_max 0 --stage_min 0 --stage_max 0 --metrics_shrink 5 \
    --input_file $tped_path --output_file_stub $comet_output_dir/out --threshold 4,4,4,3.296

echo "finished kidney_mononuclear-phagocyte in kidney"