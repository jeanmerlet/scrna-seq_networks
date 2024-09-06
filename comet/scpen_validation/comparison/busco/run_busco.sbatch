#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J pfc_glial_non_neuron_micro_output_sweep
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 12:00:00
#SBATCH -p batch
#SBATCH -N 1


source activate /lustre/orion/syb112/proj-shared/environments/andes/rwr_toolkit

BASE_DIR=/lustre/orion/syb111/proj-shared/Personal/lanemj/projects/rwr_kfold_comparisons

NETWORK_DIR_FILEPATH=$BASE_DIR/data/pfc_glial-non-neuron_micro/output/r2_ensemble/  # Base networks you wish to run cv on
FLIST_DIR_FILEPATH=$BASE_DIR/flists/r2_ensemble   # 
MULTIPLEX_DIR_FILEPATH=$BASE_DIR/network_layers/ #
GOLDSET_DIR_FILEPATH=/lustre/orion/syb112/proj-shared/personal/lanemj/pasc/pasc_permutations_abstract/human_genes/  #
CV_OUTPUT_DIR_FILEPATH=$BASE_DIR/cv_outputs #



# 1. Create the flists for every possible network layer in the $NETWORK_DIR_FILEPATH
echo Creating Flists: 
Rscript $BASE_DIR/create_layers_as_multiplexes.R -i $NETWORK_DIR_FILEPATH -o $FLIST_DIR_FILEPATH

# 2. Create multiplex Rdata file for every flist in the FLIST_DIR_FILEPATH
echo Creating Multiplexes: 
Rscript $BASE_DIR/run_rwr_make_multiplex.R -i $FLIST_DIR_FILEPATH -o $MULTIPLEX_DIR_FILEPATH

# 3. Run RWR_CV on all genes in the GOLDSET_DIR_FILEPATH (all 300 busco genes) for each Multiplex file in the MULTIPLEX_DIR_FILEPATH 
echo Run RWR CV Kfold: 
Rscript $BASE_DIR/run_rwr_kfold.R -i $MULTIPLEX_DIR_FILEPATH -g $GOLDSET_DIR_FILEPATH -o $CV_OUTPUT_DIR_FILEPATH
