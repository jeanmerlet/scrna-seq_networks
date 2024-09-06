#!/bin/bash
#SBATCH -A SYB111
#SBATCH -J setup_irf
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -o /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/logs/setup_irf.%J.out
#SBATCH -e /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/logs/setup_irf.%J.err


source activate python_andes


cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_cge-inhibitory-interneurons_id2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_cge-inhibitory-interneurons_id2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_cge-inhibitory-interneurons_lamp5_nos1

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_cge-inhibitory-interneurons_lamp5_nos1_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_cge-inhibitory-interneurons_vip

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_cge-inhibitory-interneurons_vip_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_excitatory-principal-neurons_l2-3_cux2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_excitatory-principal-neurons_l2-3_cux2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_excitatory-principal-neurons_l4_rorb

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_excitatory-principal-neurons_l4_rorb_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_excitatory-principal-neurons_l5-6_themis

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_excitatory-principal-neurons_l5-6_themis_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_excitatory-principal-neurons_l5-6_tle4

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_excitatory-principal-neurons_l5-6_tle4_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_glial-non-neuron_astro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_glial-non-neuron_astro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_glial-non-neuron_micro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_glial-non-neuron_micro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_glial-non-neuron_oligo

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_glial-non-neuron_oligo_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_glial-non-neuron_oligo_mat

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_glial-non-neuron_oligo_mat_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_glial-non-neuron_opc

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_glial-non-neuron_opc_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

echo "done"

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_glial-non-neuron_vas

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_glial-non-neuron_vas_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_mge-inhibitory-interneurons_pv

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_mge-inhibitory-interneurons_pv_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_mge-inhibitory-interneurons_pv_scube3

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_mge-inhibitory-interneurons_pv_scube3_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba10_mge-inhibitory-interneurons_sst

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba10_mge-inhibitory-interneurons_sst_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_cge-inhibitory-interneurons_id2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_cge-inhibitory-interneurons_id2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_cge-inhibitory-interneurons_lamp5_nos1

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_cge-inhibitory-interneurons_lamp5_nos1_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_cge-inhibitory-interneurons_vip

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_cge-inhibitory-interneurons_vip_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_excitatory-principal-neurons_l2-3_cux2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_excitatory-principal-neurons_l2-3_cux2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_excitatory-principal-neurons_l4_rorb

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_excitatory-principal-neurons_l4_rorb_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_excitatory-principal-neurons_l5-6_themis

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_excitatory-principal-neurons_l5-6_themis_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_excitatory-principal-neurons_l5-6_tle4

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_excitatory-principal-neurons_l5-6_tle4_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_glial-non-neuron_astro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_glial-non-neuron_astro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_glial-non-neuron_micro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_glial-non-neuron_micro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_glial-non-neuron_oligo

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_glial-non-neuron_oligo_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_glial-non-neuron_oligo_mat

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_glial-non-neuron_oligo_mat_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_glial-non-neuron_opc

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_glial-non-neuron_opc_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_glial-non-neuron_vas

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_glial-non-neuron_vas_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_mge-inhibitory-interneurons_pv

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_mge-inhibitory-interneurons_pv_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_mge-inhibitory-interneurons_pv_scube3

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_mge-inhibitory-interneurons_pv_scube3_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba8_mge-inhibitory-interneurons_sst

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba8_mge-inhibitory-interneurons_sst_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_cge-inhibitory-interneurons_id2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_cge-inhibitory-interneurons_id2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_cge-inhibitory-interneurons_lamp5_nos1

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_cge-inhibitory-interneurons_lamp5_nos1_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_cge-inhibitory-interneurons_vip

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_cge-inhibitory-interneurons_vip_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_excitatory-principal-neurons_l2-3_cux2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_excitatory-principal-neurons_l2-3_cux2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_excitatory-principal-neurons_l4_rorb

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_excitatory-principal-neurons_l4_rorb_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_excitatory-principal-neurons_l5-6_themis

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_excitatory-principal-neurons_l5-6_themis_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_excitatory-principal-neurons_l5-6_tle4

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_excitatory-principal-neurons_l5-6_tle4_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_glial-non-neuron_astro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_glial-non-neuron_astro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_glial-non-neuron_micro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_glial-non-neuron_micro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_glial-non-neuron_oligo

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_glial-non-neuron_oligo_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_glial-non-neuron_oligo_mat

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_glial-non-neuron_oligo_mat_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_glial-non-neuron_opc

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_glial-non-neuron_opc_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_glial-non-neuron_vas

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_glial-non-neuron_vas_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_mge-inhibitory-interneurons_pv

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_mge-inhibitory-interneurons_pv_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_mge-inhibitory-interneurons_pv_scube3

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_mge-inhibitory-interneurons_pv_scube3_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/ba9_mge-inhibitory-interneurons_sst

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_ba9_mge-inhibitory-interneurons_sst_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/cge-inhibitory-interneurons_id2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_cge-inhibitory-interneurons_id2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/cge-inhibitory-interneurons_lamp5_nos1

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_cge-inhibitory-interneurons_lamp5_nos1_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/cge-inhibitory-interneurons_vip

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_cge-inhibitory-interneurons_vip_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/excitatory-principal-neurons_l2-3_cux2

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_excitatory-principal-neurons_l2-3_cux2_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/excitatory-principal-neurons_l4_rorb

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_excitatory-principal-neurons_l4_rorb_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/excitatory-principal-neurons_l5-6_themis

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_excitatory-principal-neurons_l5-6_themis_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/excitatory-principal-neurons_l5-6_tle4

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_excitatory-principal-neurons_l5-6_tle4_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/glial-non-neuron_astro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_glial-non-neuron_astro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/glial-non-neuron_micro

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_glial-non-neuron_micro_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/glial-non-neuron_oligo

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_glial-non-neuron_oligo_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/glial-non-neuron_oligo_mat

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_glial-non-neuron_oligo_mat_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/glial-non-neuron_opc

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_glial-non-neuron_opc_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/glial-non-neuron_vas

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_glial-non-neuron_vas_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/mge-inhibitory-interneurons_pv

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_mge-inhibitory-interneurons_pv_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/mge-inhibitory-interneurons_pv_scube3

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_mge-inhibitory-interneurons_pv_scube3_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

cd /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/code/irf-loop/irf-runs/brain/pfc/mge-inhibitory-interneurons_sst

python /gpfs/alpine/syb105/proj-shared/Projects/iRF/iRF_LOOP_SetUp.py /gpfs/alpine/syb105/proj-shared/Projects/scRNA-seq/data/human/brain/healthy/lister_lab/irf-loop_mtx/pfc_mge-inhibitory-interneurons_sst_irf-loop-mtx.tsv --TotalNodes 50 --NodesPer 10 --RunTime 2 --Account SYB111

