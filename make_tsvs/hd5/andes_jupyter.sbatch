#!/bin/bash

#SBATCH -A SYB111
#SBATCH -t 48:00:00
#SBATCH -J notebook
#SBATCH -o ./logs/notebook.%J.err
#SBATCH -e ./logs/notebook.%J.out
#SBATCH -N 1
#SBATCH -p gpu

source activate cytograph
#module load python
#jupyter nbextension enable --py widgetsnbextension

jupyter lab --no-browser --port=9009
