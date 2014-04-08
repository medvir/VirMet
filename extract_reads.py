#!/opt/python2.7/bin/python2.7
import sys
import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

ids_file = sys.argv[1]
ids = set([l.strip() for l in open(ids_file)])

fastq_file = sys.argv[2]
if fastq_file.endswith('.gz'):
    fq_handle = gzip.open(fastq_file)
else:
    fq_handle = open(fastq_file)

output_handle = open('matching.fastq', 'w')

c = 0
# Using FastqGeneralIterator allows fast performance
for title, seq, qual in FastqGeneralIterator(fq_handle):
    if title.split()[0] in ids:
        c += 1
        output_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
        if c % 10000 == 0:
            print >> sys.stderr, 'written %d clean reads' % c

output_handle.close()