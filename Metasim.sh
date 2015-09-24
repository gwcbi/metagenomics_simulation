#!/bin/bash

# REMEMBER TO LAUNCH IT IN THE METASIM EXECUTABLE FOLDER!!!!

#SBATCH -J Metasim
#SBATCH -o Metasim_%A_%a.out
#SBATCH -e Metasim_%A_%a.err
#SBATCH -t 2-00:00:00
#SBATCH --partition=2tb,defq,short,debug,gpu
#SBATCH --nodes=1

# The simulation is setup as follows:
# - The median insert length of PE reads is 550nt;
# - The probability of generating PE reads is 1 (ie, all simulated reads will be PE)
# - Reads are simulated error-free
# - 16 processors are used
#
#./MetaSim cmd -r5000000 -f550 -t100 -m --empirical-pe-probability 1.0 -g fakeError300.mconf -2 fakeError300.mconf --threads 16 ../metasim_data/reference_dataset0001_v20Evenness.mprf

working_dir=$1

cat ${working_dir}/metasim_simulation_table_1 | while read LINE
do
	datasetID=$(echo ${LINE} | awk '{print $1}')
	n_reads=$(echo ${LINE} | awk '{print $2}')
	read_length=$(echo ${LINE} | awk '{print $3}')
	abundance_file=$(echo ${LINE} | awk '{print $4}')

	echo "option row is" ${LINE}
	echo "datasetID is" ${datasetID}
	echo "number of reads is" ${n_reads}
	echo "read length is" ${read_length}

	mkdir -p ${working_dir}/${datasetID}
	rm ${working_dir}/${datasetID}/*

	./MetaSim cmd -r${n_reads} -f550 -t100 -m --empirical-pe-probability 1.0 -g /home/domenico_simone/metasim/fakeError${read_length}.mconf -2 /home/domenico_simone/metasim/fakeError${read_length}.mconf --threads 16 -z -d ${working_dir}/${datasetID} ${working_dir}/${abundance_file}

	echo "output folder is" ${datasetID}

	outfile=$(echo ${abundance_file} | sed 's/\.mprf/\-Empirical.fna.gz/g')
	mv ${working_dir}/${datasetID}/${outfile} ${working_dir}/${datasetID}/${datasetID}-reads.fa.gz
done

