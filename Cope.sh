#!/bin/bash

# MEMO: 
# - delete SBATCH args;
# - set up working directory;
# - export cope PATH

#SBATCH -J Cope
#SBATCH -o Cope_%A_%a.out
#SBATCH -e Cope_%A_%a.err
#SBATCH -t 4:00:00
#SBATCH --partition=2tb,defq,short,gpu,debug
#SBATCH --nodes=1

for dir in $(ls -d dataset_000*); do cd ${dir}; cope -a ${dir}-reads.1.fa -b ${dir}-reads.2.fa -o ${dir}-reads.merged.fa -2 ${dir}-reads.R1.fa -3 ${dir}-reads.R2.fa -m 0 >cope.log 2>cope.error; cd ..; done
