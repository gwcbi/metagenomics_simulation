#!/bin/bash

#SBATCH -J cope_150610
#SBATCH -o cope_150610_%A_%a.out
#SBATCH -e cope_150610_%A_%a.err
#SBATCH -t 4:00:00
#SBATCH --partition=2tb,defq,short,gpu,debug
#SBATCH --nodes=1
#MBATCH --array=5
#SBATCH --mail-type=ALL
#SBATCH --mail-user="domenico_simone@gwu.edu"

cd /home/domenico_simone/metasim_data/150730_3

for dir in $(ls -d dataset_000*); do cd ${dir}; cope -a ${dir}-reads.1.fa -b ${dir}-reads.2.fa -o ${dir}-reads.merged.fa -2 ${dir}-reads.R1.fa -3 ${dir}-reads.R2.fa -m 0 >cope.log 2>cope.error; cd ..; done
