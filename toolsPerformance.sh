#!/usr/bin/env bash

#SBATCH -J performance
#SBATCH -o performance_%A_%a.out
#SBATCH -e performance_%A_%a.err
#SBATCH -t 4:00:00
#SBATCH --partition=2tb,defq,short
#SBATCH --nodes=1

outfile="PathoScope_performance.tsv"
dates_file="dates"

for i in $(cat ${dates_file}); do
	for folder in $(ls -d ${i}/dataset_0*/theta0); do
		echo ${folder}
		echo $(echo ${folder} | awk 'BEGIN{FS="/";OFS="\t"}{print $1, $2, $3}') $(toolsPerformance.py pathoscope ${folder} | grep -A1 Benchmark_species | tail -n1) >> ${outfile}
	done
done

#for folder in $(ls -d 150611/dataset_0*/theta0); do
##echo ${folder}
##	h=$(echo ${folder})
#echo $(echo ${folder} | awk 'BEGIN{FS="/";OFS="\t"}{print $1, $2, $3}')
#done
