#!/usr/bin/env python

# USAGE
# python performance.py <tool> <working_dir>, e.g.
# python performance.py pathoscope|pathoscope_reads dataset_0001/theta0 (PathoScope)
# python performance.py metamix dataset_0001/metaMix path/to/blastdb (metaMix)

# the script expects to find this data structure:
#
# PathoScope:
# dataset_000?
#	|- dataset_000?.-reads.fa.gz
#	|- theta?
#		|- dataset_000?.theta0-sam-report.tsv
#		|- updated_dataset_000?.sam.gz
#
# metaMix:
# dataset_000?
#	|- dataset_000?.-reads.fa.gz
#	|- metaMix
#		|- res
#			|- presentSpecies_assignedReads.tsv


# sam file has headers with @, skip them
# analyze rows like
# r2.1    83      gi|197261329|ref|NC_011165.1|   3939    255     75M     =       3805    -209    GGATGGCTGAGCTTGGGCTCTACGGGAACGGCCGGAGCGCTAGCGCTGCAACCGGCGTCTCCTGGGCCGTAGCCG     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     AS:i:150        XS:i:150
#         XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:150        YT:Z:CP
# r2.2    163     gi|197261329|ref|NC_011165.1|   3805    255     75M     =       3939    209     AGTAAAGCATTTTCTTTCATTATTTCGGAAGAGCCTGGAAATAGTTCCAGATCCAATCGCCTGCGGCCAGCAAGA     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     AS:i:150        XS:i:132
#         XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:150        YT:Z:CP
# r3.1    83      gi|197261329|ref|NC_011165.1|   8322    138     75M     =       8169    -228    AATCCCTATGTAGTTCCGACCATGCTGCAGGACTGGTATAACTCCCAAGGATTCATCGGATACCAAGCTTGCGCC     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     AS:i:150        XS:i:150
#         XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:150        YT:Z:CP
# r3.2    163     gi|197261329|ref|NC_011165.1|   8169    138     75M     =       8322    228     CTTGGGAAGATTCGCGGCTGGAACGTCGAGCCGGAGAAAGCTCCGGTCATCCGTAGCGTGAAGGATTTTCTGGAG     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII     AS:i:150        XS:i:150


# original mapping data of simulated reads have to be taken
# from the original simulation 
import sys, gzip, os, shlex, subprocess

def proc_sam_row(l): # line is not splitted
	try:
		l = l.strip().split()
		readID = l[0]
		mapping_gi = l[2].split('|')[1]
	except:
		print l, "raised an error"
	return readID, mapping_gi

def get_read_origin(l):
	# read heading as output by metasim is
	# >r1.1 |SOURCES={GI=703773128,bw,6419-6494}|ERRORS={}|SOURCE_1="Zaire ebolavirus isolate H.sapiens-tc/GIN/14/WPG-C15" (1894b8c11d315380698eec5901681ffb4fe3219c)\n'
	readID = l.split()[0][1:].strip()
	mapping_gi = l.split('GI=')[1].split(',')[0]
	return readID, mapping_gi

def genome_stats(readID, mapping_gi):
	c = 0 # number of reads from genome gi

def comparisons(benchmark, test, total_genome_list, read_analysis=False): # dictionaries {gi : reads}
	# number of correct genomes between test and benchmark set
	print "Benchmark: ", set(benchmark.keys())
	print "Test: ", set(test.keys())
	n_correct_gi = len(set(benchmark.keys()) & set(test.keys())) # tp
	n_real_gi = len(set(benchmark.keys()))
	n_test_gi = len(set(test.keys()))
	fp_gi = len(set(test.keys())) - n_correct_gi
	fn_gi = len(set(benchmark.keys())) - n_correct_gi
	tn_gi = len(set(total_genome_list) - (set(benchmark.keys()) | set(test.keys())))
	sn_gi = float(n_correct_gi)/(n_correct_gi + fn_gi)
	sp_gi = float(tn_gi)/(tn_gi + fp_gi)
	# ji as tp/(tp + fn + fp)
	ji_gi = float(n_correct_gi)/(n_correct_gi + fn_gi + fp_gi)
	print "\nStatistics on number of identified species"
	print "\t".join(["Benchmark_species", "Identified_species", "Correctly_identified_species", "Sensitivity", "Specificity", "Jaccard_index"])
	print "\t".join([str(i) for i in [len(benchmark.keys()), len(test.keys()), n_correct_gi, str(sn_gi)[:4], str(sp_gi)[:4], str(ji_gi)[:4]]])
#	print "Total species in benchmark: {0}".format(str(len(benchmark.keys())))
#	print "Total species identified: {0}".format(str(len(test.keys())))
#	print "Total species correctly identified: {0}".format(str(n_correct_gi))
#	print "Sensitivity: {0}".format(str(sn_gi)[:4])
#	print "Specificity: {0}".format(str(sp_gi)[:4])
#	print "Jaccard Index: {0}".format(str(ji_gi)[:4])

	if read_analysis == True:
		# calculate total reads
		total_reads = sum(len(benchmark[gi]) for gi in benchmark.keys())
		# number of correct sequences between test and benchmark set
		print "\nStatistics on assigned reads based on %d genomes (out of %d identified)" % (n_correct_gi, len(test.keys()))
		print "Total reads: %d" % total_reads
		# table header
		print "\t".join(["Genome_gi", "Simulated_reads", "Assigned_reads", "Correctly_assigned_reads", "Sensitivity", "Specificity", "Jaccard_index"])
		non_test_gi = []
		for gi in benchmark.keys():
			if gi in test.keys():
				n_correct_reads = len(set(benchmark[gi]) & set(test[gi])) # tp
				n_real_reads = len(set(benchmark[gi]))
				n_test_reads = len(set(test[gi]))
				fp_reads = len(set(test[gi])) - n_correct_reads
				fn_reads = len(set(benchmark[gi])) - n_correct_reads
				tn_reads = total_reads - len(set(benchmark[gi]) | set(test[gi]))
				sp_reads = float(tn_reads)/(tn_reads + fp_reads)
				sn_reads = float(n_correct_reads)/(n_correct_reads + fn_reads) # specificity as n_correct/n_test
				ji_reads = float(n_correct_reads)/(n_correct_reads + fn_reads + fp_reads)
				print "\t".join([str(i) for i in [gi, n_real_reads, n_test_reads, n_correct_reads, str(sn_reads)[:4], str(sp_reads)[:4], str(ji_reads)[:4]]])
			else:
				non_test_gi.append(gi)
		print "\nGenomes not identified:", ",".join(non_test_gi)
	else:
		print "\nGenomes not identified:", ",".join(set(benchmark.keys()) - set(test.keys()))

def metaMix_res_processing(over_threshold_file_in, path_to_blastdb):
	"""
		Input:
		- result file: presentSpecies_assignedReads.tsv
		- path to blastdb
		
		Output:
		- total_genome_list: a list of all genomes [eg in bowtie index (PathoScope), blastdb (metaMix)]
		- d_test: dictionary { gi : [read1, read2 ...]} metaMix DOES NOT provide output mapping reads but it's for compatibility with PathoScope result processing
	"""
	total_genome_list = []
	d_test = {}
	over_threshold_file = open(over_threshold_file_in, 'r')
	r = over_threshold_file.readline()
	while r:
		if r[1].isdigit():
			d_test[shlex.split(r)[1]] = []
		r = over_threshold_file.readline()
	cmd = ['blastdbcmd', '-db', path_to_blastdb, '-outfmt', '%g', '-entry', 'all']
	blastdbcmd_query = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	o, e = blastdbcmd_query.communicate()
	total_genome_list = o.split()
	return total_genome_list, d_test

def pathoscope_res_processing(over_threshold_file_in, test_file_in, read_analysis=False):
	"""
		Input:
		- result file, eg dataset_0001.theta0-sam-report.tsv
		- sam file, eg updated_dataset_0001.sam.gz
		
		Output:
		- total_genome_list: a list of all genomes [eg in bowtie index (PathoScope), blastdb (metaMix)]
		- d_test: dictionary { gi : [read1, read2 ...]}
	"""
	total_genome_list = [] # count total genomes
	try:
		test_file = open(test_file_in, 'r')
	except:
		test_file = gzip.open(test_file_in+".gz")
	d_test = {}
	over_threshold_file = open(over_threshold_file_in, 'r')
	r = over_threshold_file.readline()
	while r:
		if r.startswith('gi'):
			d_test[r.split('|')[1]] = []
		r = over_threshold_file.readline()
	if read_analysis == True:
		for k in test_file:
			if k[0] != '@':
				readID, mapping_gi = proc_sam_row(k)
				if mapping_gi in d_test.keys():
					d_test[mapping_gi].append(readID)
				#else:
				#	d_test[mapping_gi] = [readID]
			elif k.startswith('@SQ'):
				total_genome_list.append(k.split('|')[1])
	else:
		for k in test_file:
			if k.startswith('@SQ'):
				total_genome_list.append(k.split('|')[1])
	return total_genome_list, d_test

tool = sys.argv[1]
if tool.lower() == "pathoscope":
	analysis = "PathoScope"
elif tool.lower() == "pathoscope_reads":
	analysis = "PathoScope_Reads"
elif tool.lower() == "metamix":
	analysis = "metaMix"
else:
	sys.exit('No valid tool specified.')

# /path/to/metasim/data/150611/dataset_0001/theta0 (Pathoscope)
# /path/to/metasim/data/150611/dataset_0001/metaMix/res
working_dir = os.path.abspath(sys.argv[2])

working_dir_supfolder, working_dir_folder = os.path.split(working_dir) # /path/to/metasim/data/150611/dataset_0001/
dataset_name = os.path.split(working_dir_supfolder)[1] # dataset_0001
benchmark_file_in = os.path.join(working_dir_supfolder, dataset_name+'-reads.fa.gz') # dataset_0001

# process benchmark, ie the compressed fasta read file. Common to all tools performance evaluation
try:
	benchmark_file = gzip.open(benchmark_file_in)
except:
	benchmark_file = open(benchmark_file_in, 'r')
d_bench = {}
for l in benchmark_file:
	if l[0] == '>':
		readID, mapping_gi = get_read_origin(l)
		if mapping_gi in d_bench.keys():
			d_bench[mapping_gi].append(readID)
		else:
			d_bench[mapping_gi] = [readID]

# run the analysis
if analysis == "PathoScope":
	test_file_in = os.path.join(working_dir, "updated_"+dataset_name+".sam") # updated_dataset_0001.sam; if it doesn't exist, will be updated_dataset_0001.sam.gz
	over_threshold_file_in = os.path.join(working_dir, dataset_name+"."+working_dir_folder+"-sam-report.tsv")
	total_genome_list, d_test = pathoscope_res_processing(over_threshold_file_in, test_file_in)
	comparisons(d_bench, d_test, total_genome_list)
elif analysis == "PathoScope_Reads":
	test_file_in = os.path.join(working_dir, "updated_"+dataset_name+".sam") # updated_dataset_0001.sam; if it doesn't exist, will be updated_dataset_0001.sam.gz
	over_threshold_file_in = os.path.join(working_dir, dataset_name+"."+working_dir_folder+"-sam-report.tsv")
	total_genome_list, d_test = pathoscope_res_processing(over_threshold_file_in, test_file_in, read_analysis=True)
	comparisons(d_bench, d_test, total_genome_list, read_analysis=True)
elif analysis == "metaMix":
	path_to_blastdb = sys.argv[3]
	over_threshold_file_in = os.path.join(working_dir, 'res', 'presentSpecies_assignedReads.tsv')
	total_genome_list, d_test = metaMix_res_processing(over_threshold_file_in, path_to_blastdb)
	comparisons(d_bench, d_test, total_genome_list)



