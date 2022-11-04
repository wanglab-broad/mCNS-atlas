#!/bin/bash

# Request memory
#$ -l h_vmem=200G

# Request runtime
#$ -l h_rt=05:00:00

# Specify output file destination
#$ -o /stanley/WangLab/kamal/outputs/integration/spatial/

# Join output and error files
#$ -j y


# Set up Python environment
source /broad/software/scripts/useuse
reuse Anaconda3
source activate /stanley/WangLab/envs/neighbors

# Args
num_nbrs=$(($SGE_TASK_ID*10))

# Run script
python /stanley/WangLab/kamal/code/integration/spatial/level2_k_titration.py --num_nbrs $num_nbrs