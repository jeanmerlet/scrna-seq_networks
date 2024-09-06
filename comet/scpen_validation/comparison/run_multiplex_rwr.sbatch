#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J rwr
#SBATCH -o logs/rwr.%j.out
#SBATCH -e logs/rwr.%j.err
#SBATCH -t 12:00:00
#SBATCH -p batch
#SBATCH -N 1

source activate /lustre/orion/syb112/proj-shared/environments/andes/rwr_toolkit

base_dir="/lustre/orion/syb111/proj-shared/Personal/lanemj/projects/rwr_kfold_comparisons"
networks_dir="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_out/pfc_glial-non-neuron_micro/networks/networks"
flist_dir="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_out/pfc_glial-non-neuron_micro/networks/flist"
multiplex_dir="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_out/pfc_glial-non-neuron_micro/networks/multiplex"
goldset_dir="/lustre/orion/syb112/proj-shared/personal/lanemj/pasc/pasc_permutations_abstract/human_genes/"
rwr_output_dir="/gpfs/alpine2/syb112/proj-shared/Projects/scrna-seq/data/human/brain/healthy/lister_lab/comet_out/pfc_glial-non-neuron_micro/networks/rwr"


echo Creating Flists: 
#Rscript $base_dir/create_layers_as_multiplexes.R -i $networks_dir -o $flist_dir

echo Creating Multiplexes: 
Rscript $base_dir/run_rwr_make_multiplex.R -i $flist_dir -o $multiplex_dir

echo Run RWR CV Kfold: 
Rscript $base_dir/run_rwr_kfold.R -i $multiplex_dir -g $goldset_dir -o $rwr_output_dir
