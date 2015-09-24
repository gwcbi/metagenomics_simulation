# USAGE
# python SplitMetasimPairs.py <folder> (with name consistent with compressed fasta inside)
# eg python SplitMetasimPairs.py dataset_0001
# will analyze the file dataset_0001/dataset_0001-reads.fa.gz

import gzip, sys
from Bio import SeqIO

input_folder = sys.argv[1]
input_file = '/'.join([sys.argv[1], sys.argv[1]+"-reads.fa.gz"])

h = gzip.open(input_file, 'r')
d1 = open(input_file.replace("reads.fa.gz", "reads.1.fa"), 'w')
#d1 = open('d_R1.fa', 'w')
d2 = open(input_file.replace("reads.fa.gz", "reads.2.fa"), 'w')
for i in SeqIO.parse(h, 'fasta'):
	if i.id[-2:] == '.1':
		d1.write('>%s\n%s\n' % (i.id, str(i.seq)))
	elif i.id[-2:] == '.2':
		d2.write('>%s\n%s\n' % (i.id, str(i.seq)))

d1.close()
d2.close()
