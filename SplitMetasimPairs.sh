#!/bin/bash
#SBATCH -J SplitMetasimPairs
#SBATCH -o SplitMetasimPairs_%A_%a.out
#SBATCH -e SplitMetasimPairs_%A_%a.err
#SBATCH -t 4:00:00
#SBATCH --partition=2tb,defq,short,debug
#SBATCH --array=1-9
#SBATCH --nodes=1

input_name=$(sed -n "$SLURM_ARRAY_TASK_ID"p dataset_file)

python SplitMetasimPairs.py ${input_name}
